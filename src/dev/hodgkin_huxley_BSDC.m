%HODGKIN_HUXLEY_BSDC
%   V = HODGKIN_HUXLEY_DYNAMICCLAMP(T,I_EXT) simulates the membrane voltage 
%   V of the Hodgkin-Huxley neuron model in response to an external current
%   I_EXT over times T. The lengths of T and I_EXT must be the same. 
%   
%   V = HODGKIN_HUXLEY_DYNAMICCLAMP(...,'MEMCON',@MEM_CON) injects a
%   memconductance into the Hodgkin-Huxley neuron model, and simulates the 
%   membrane voltage V. The memconductance is described by a set of 
%   ordinary differential equations, which are implemented in the function 
%   handle MEM_CON. User should also use 
%   HODGKIN_HUXLEY_DYNAMICCLAMP(...,'MEMINITSTATE',STATE)
%   to initialize the value of the internal state of the memconductance. 
%
%   HODGKIN_HUXLEY_DYNAMICCLAMP(...,'MAXCON',MAXCON,'MINCON',MINCON) 
%   specifies the initial value for the upper and lower bound of the maximum
%   conductance value.
%
%   HODGKIN_HUXLEY_DYNAMICCLAMP(...,'MAXVOLT',MAXVOLT,'MINVOLT',MINVOLT) 
%   sets the boundary values for determining whether the system is brokendown.
% 
%   HODGKIN_HUXLEY_DYNAMICCLAMP(...,'RESETPERIOD',RESETPERIODTIME) sets 
%   RESETRERIOD, the duration of period of time, in which the system enters the 
%   reset mode. In reset mode the value of conductance of the injected devise is 
%   set to zero. The system should return back to the steady state of zero-input 
%   response. If not, RESETPERIOD needs to be extended.
%
%   HODGKIN_HUXLEY_DYNAMICCLAMP(...,'WAITPERIOD',WAITPERIOD) sets WAITPERIOD, 
%   the maximal duration of time that the system is not brokendown. This 
%   parameter is used to judge that the system is in a stable mode.
%
%   V = HODGKIN_HUXLEY_DYNAMICCLAMP(...,'HHN',@hhn) replaces the function 
%   handle for Hodgkin-Huxley equations with the @HHN, and then simulates 
%   the membrane voltage V.
%
%   HODGKIN_HUXLEY_DYNAMICCLAMP(...,'VERBOSE',TRUE) sets the program to verbose
%   mode. Internal state of the program will be printed at command line. 
%   'VERBOSE' is set to false by default.
%
%   Authors: Chung-Heng Yeh <chyeh@ee.columbia.edu>
%
%   Copyright 2012-2014 Aurel A. Lazar, Yiyin Zhou, and Chung-Heng Yeh

function [Vout Mout ctrl_signal] = hodgkin_huxley_BSDC(t, I_ext, varargin)
    
    % Handle the optional input parameters.
    % =====================================================================
    p = inputParser;
    % Specify the Hodgkin-Huxley ODEs using a function handle.
    addParamValue(p,'HHN',@default_hhn,@(x) isa(x,'function_handle'));
    % Turn 'On' or 'Off' the sodium, potassium, and leaky channel.
    addParamValue(p,'Sodium',   'on', @(x) any(validatestring(x,{'On','Off'})));
    addParamValue(p,'Potassium','on', @(x) any(validatestring(x,{'On','Off'})));
    addParamValue(p,'Leaky',    'on', @(x) any(validatestring(x,{'On','Off'})));
    
    % Specify the memconductance ODEs using a function handle, and the 
    % memconductance initial state.
    addParamValue(p,'MemInitState', [0 0 0],@isnumeric);
    null_mem_con = @(v,x) [0 x.*[0 0 0]];  
    addParamValue(p,'MemCon', null_mem_con,@(x) isa(x,'function_handle'));
    
    % Specify the maximum conductance value
    addParamValue(p,'MaxCon', null_mem_con, @isnumeric);
    % Specify the minimum conductance value
    addParamValue(p,'MinCon', Inf, @isnumeric);
    % Sspecify the maximum conductance update step
    %addParamValue(p,'MaxConUpdate', Inf, @isnumeric);

    % Specify the upper bound of the membrane voltage
    addParamValue(p,'MaxVolt',   100, @isnumeric);
    % Specify the lower bound of the membrane voltage
    addParamValue(p,'MinVolt', -100, @isnumeric);
    % Specify the reset period
    addParamValue(p,'resetPeriod', .05, @isnumeric);
    % Specify the waiting period
    addParamValue(p,'waitPeriod', .2, @isnumeric);
    
    % Set/Unset Verbose mode
    addParamValue(p,'Verbose',false,@islogical);
    
    p.KeepUnmatched = true;
    parse(p,varargin{:});

    % Handle the unexpected input parameter.
    UnmatchedParam = fieldnames(p.Unmatched);
    if ~isempty(UnmatchedParam)
        error(['"',UnmatchedParam{1},'" is not a valid parameter.']);
    end
    % =====================================================================
    % Set display mode
    if p.Results.Verbose
        mydisp = @(varargin) fprintf(varargin{:});
    else
        mydisp = @(varargin) true;
    end
    % Use boolean flag to indicate the on/off status of each channel.
    NaKL_flag = strcmpi('on', ...
                {p.Results.Sodium, p.Results.Potassium, p.Results.Leaky});

    % Assume that the time is given in seconds, convert time to millisecond.
    dt = 1000*(t(2)-t(1));
    
    % Initialize HHN states V, n, m, and h.
    nmh_state = [0 0 1];
    v    = -64.9964;
    V3   = zeros(1,3);
    Vout = zeros(numel(t),4);
    spk_status = zeros(1,4);
    % Initialize the Hodgkin-Huxley neuron function handel. 
    hhn = p.Results.HHN;
    
    % Initialize Memconductance internal states and function handel. 
    % Notice that the output of mem_ode is [Im  d_mem_state]. Im is the 
    % term in V's ODE caused by the memristor.
    mem_state   = p.Results.MemInitState;
    mem_con     = p.Results.MemCon;
    mem_max_con = p.Results.MaxCon*2; % Scaled by 2 here, but would be scaled down
    mem_min_con = p.Results.MinCon;
    %mem_max_up  = p.Results.MaxConUpdate;
    
    % Convert reset period to number of steps
    reset_step  = round((p.Results.resetPeriod)*1000/dt );
    reset_index = 1+reset_step;
    reset_flag  = true;
    % Convert waiting period to number of steps
    wait_step   = round((p.Results.waitPeriod)*1000/dt );
    wait_index  = reset_index + wait_step;
    % Initialize mem conductance value
    MemConVal = zeros(size(mem_max_con));
    Mout      = zeros(numel(t),numel(MemConVal));
    
    % Initialize contral signal vector; column-1 is for reset signal, and 
    % column-2 is for breakdown time.
    ctrl_signal = zeros(numel(t),2);
    
    % Use forward Euler method to solve ODE.
    for i = 1:numel(I_ext)
        
        % Leave Reset mode
        if i == reset_index
            MemConVal  = 0.5*( mem_max_con + mem_min_con );
            reset_flag = false;
            mydisp('Conductance is set to: %s...',num2str(MemConVal));
        end
        ctrl_signal(i,1) = reset_flag;
        
        Mout(i,:) = MemConVal; 
        % Compute the HHN ODEs. Notice the first entry of the output is 
        % used in dv/dt equitoan.
        temp  = hhn(v, nmh_state, NaKL_flag);
        I_hhn = temp(1);
        d_nmh = temp(2:end);
        % Compute the Memconductance ODEs. Notice the first entry of the output 
        % is used in dv/dt equation.
        temp  = mem_con(v, mem_state, MemConVal);
        I_mem = temp(1);
        d_mem = temp(2:end);
        % Use the forward Euler method to integrate.
        mem_state  = mem_state + dt * d_mem;
        nmh_state  = nmh_state + dt * d_nmh;
        dv = (I_ext(i)+I_hhn+I_mem);
        v  = v + dt * dv;
        V3 = [V3(2:3) v];
        Vout(i,:) = [v nmh_state];
        %Vout(i,1) = dv;
        if ~reset_flag
            % Membrane voltage tends to be unbounded, indicating system is brokendown
            if v > p.Results.MaxVolt || v < p.Results.MinVolt
                mydisp('system breaks down.\n');
                mem_max_con = MemConVal;
                reset_mode();
                ctrl_signal(i,2) = 1;
            end
            % Membrane voltage reaches steady state
            if norm( [dv d_nmh] ) < 1e-10
                mydisp('reaches the steady state.\n');
                mem_min_con = MemConVal;
                reset_mode();
            end
            % Detecting Spike
            if V3(2) > V3(1) && V3(2) > V3(3)
                if norm( (spk_status-Vout(i-1,:))./spk_status ) < 5e-3
                    mydisp('operating on the limit cycle.\n');
                    mem_min_con = MemConVal;
                    reset_mode();
                else
                    spk_status = Vout(i-1,:);
                end
            end
        end
    end
    mydisp('End of dynamic clamping\n');
    function reset_mode
        %mem_max_con = mem_max_con + mem_max_up;
        MemConVal   = zeros(size(MemConVal));
        reset_index = i+reset_step;
        wait_index  = reset_index + wait_step;
        reset_flag  = true;
        spk_status  = ones(size(spk_status))*-10;
    end
end
% Define Default Hodgkin-Huxley Nueron ODEs. This function handle is 
% written by Yiyin Zhou.
function dx = default_hhn(v,x,NaKL_flag)
    g = [120 36 0.3].*NaKL_flag;
    a = exp(-(v+55)/10)-1;
    if a == 0
        dn = (1-x(1)) * 0.1 - x(1) * (0.125*exp(-(v+65)/80));
    else
        dn = (1-x(1)) * (-0.01*(v+55)/a) - x(1) * (0.125*exp(-(v+65)/80));
    end

    a = exp(-(v+40)/10)-1;
    if a == 0
        dm = (1-x(2)) - x(2) * (4*exp(-(v+65)/18));
    else
        dm = (1-x(2)) * (-0.1*(v+40)/a) - x(2) * (4*exp(-(v+65)/18));
    end
    dh = (1-x(3)) * (0.07*exp(-(v+65)/20)) - x(3) / (exp(-(v+35)/10)+1);
    I_nmh = - g(1)*x(2).^3*x(3)*(v-50) - g(2)*x(1).^4*(v+77) - g(3)*(v+54.387);
    dx    = [I_nmh dn dm dh].*[1 NaKL_flag(2) NaKL_flag(1) NaKL_flag(1)];    
end