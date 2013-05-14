%HODGKIN_HUXLEY_dynamicClamp 
%   V = HODGKIN_HUXLEY_DYNAMICCLAMP(T,I_EXT) simulates the membrane voltage 
%   V of the Hodgkin-Huxley neuron model in response to an external current
%   I_EXT over times T. The lengths of T and I_EXT must be the same. 
% 
%   V = HODGKIN_HUXLEY_DYNAMICCLAMP(...,Channel,'On'/'Off',...) turns ON or
%   OFF the CHANNEL in the Hodgkin-Huxley neuron model, and simulates the 
%   membrane voltage V in response to an external current I_EXT over times 
%   T. The CHANNEL can be 'SODIUM', 'POTASSIUM', or 'LEAKY'. When no 
%   CHANNEL is specified, every channel is turned on. User can specify the 
%   status of mutiple channels simultaneously, for example,
%   HODGKIN_HUXLEY_DYNAMICCLAMP(...,'LEAKY','OFF','POTASSIUM','OFF').
%   
%   V = HODGKIN_HUXLEY_DYNAMICCLAMP(...,'MEMCON',@MEM_CON) injects a
%   memconductance into the Hodgkin-Huxley neuron model, and simulates the 
%   membrane voltage V. The memconductance is described by a set of 
%   ordinary differential equations, which are implemented in the function 
%   handle MEM_CON. User should also use 
%   HODGKIN_HUXLEY_DYNAMICCLAMP(...,'MEMINITSTATE',STATE)
%   to initialize the value of the internal state of the memconductance. 
%
%   V = HODGKIN_HUXLEY_DYNAMICCLAMP(...,'MAXCON',maxcon) specifies the 
%   maximum conductance of the injected mem-conductance.
%
%   V = HODGKIN_HUXLEY_DYNAMICCLAMP(...,'ResetTime',TIME) sets the TIME 
%   at which the maximum values of mem-conductance is applied. Before TIME, 
%   the maximum conductance of the mem-conductance is initialized to a zero 
%   vector. If 'ResetTime' is not specified, TIME is set to the beginning 
%   of the simulation.
%
%   V = HODGKIN_HUXLEY_DYNAMICCLAMP(...,'HHN',@hhn) replaces the function 
%   handle for Hodgkin-Huxley equations with the @HHN, and then simulates 
%   the membrane voltage V.
%
%   This file is based on the code written by Yiyin Zhou.    
%
%   Authors: Chung-Heng Yeh
%
%   Copyright 2012-2014 Aurel A. Lazar, Yiyin Zhou, and Chung-Heng Yeh

function [Vout] = hodgkin_huxley_dynamicClamp(t, I_ext, varargin)
    
    % Handle the optional input parameters. You can simply jump to line 59.
    % =====================================================================
    p = inputParser;
    % Specify the Hodgkin-Huxley ODEs using a function handle.
    addParamValue(p,'HHN',@default_hhn,@(x) isa(x,'function_handle'));
    % Specify the memconductance ODEs using a function handle, and the 
    % memconductance initial state.
    addParamValue(p,'MemInitState', [0 0 0],@isnumeric);
    null_mem_con = @(v,x) [0 x.*[0 0 0]];  
    addParamValue(p,'MemCon', null_mem_con,@(x) isa(x,'function_handle'));
    % Specify the time to reset the maximum value of the mem-conductance
    addParamValue(p,'resetTime', t(1),@isnumeric);
    % Specify the maximum conductance value
    addParamValue(p,'MaxCon', null_mem_con,@isnumeric);
    % Turn 'On' or 'Off' the sodium, potassium, and leaky channel.
    addParamValue(p,'Sodium',   'on', @(x) any(validatestring(x,{'On','Off'})));
    addParamValue(p,'Potassium','on', @(x) any(validatestring(x,{'On','Off'})));
    addParamValue(p,'Leaky',    'on', @(x) any(validatestring(x,{'On','Off'})));

    p.KeepUnmatched = true;
    parse(p,varargin{:});

    % Handle the unexpected input parameter.
    UnmatchedParam = fieldnames(p.Unmatched);
    if ~isempty(UnmatchedParam)
        error(['"',UnmatchedParam{1},'" is not a valid parameter.']);
    end
    % =====================================================================
    
    % Use boolean flag to indicate the on/off status of each channel.
    NaKL_flag = strcmpi('on', ...
                {p.Results.Sodium, p.Results.Potassium, p.Results.Leaky});

    % Assume that the time is given in seconds, convert time to millisecond.
    dt = 1000*(t(2)-t(1));
    
    % Initialize HHN states V, n, m, and h.
    nmh_state = [0 0 1];
    v = -64.9964;
    Vout = zeros(numel(t),4);
    
    % Initialize the Hodgkin-Huxley neuron function handel. 
    hhn = p.Results.HHN;
    
    % Initialize Memconductance internal states and function handel. 
    % Notice that the output of mem_ode is [Im  d_mem_state]. Im is the 
    % term in V's ODE caused by the memristor.
    mem_state   = p.Results.MemInitState;
    mem_con     = p.Results.MemCon;
    
    % Convert reset time to discrete index
    reset_index = floor((p.Results.resetTime-t(1))*1000/dt ) + 1;
    % Initialize mem conductance value
    MemConVal = [0 0 0];
    % Use forward Euler method to solve ODE.
    for i = 1:numel(I_ext)
        if i == reset_index
            MemConVal = p.Results.MaxCon;
        end;
            
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
        v = v + dt * (I_ext(i)+I_hhn+I_mem);
        Vout(i,:) = [v nmh_state];
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