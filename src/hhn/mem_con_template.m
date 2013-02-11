%MEM_CON_TEMPLATE(v, x, t, varargin)
%   dx = EM_CON_TEMPLATE(v, x, varargin) computes the derivatives of state 
%   variables X based on a set ordinary differential equations implemented 
%   inside this function at time T. Each ODE is a function of state 
%   variables X and the voltage V. The maximum conductance of each channel 
%   is given in G. Additional arguments can be provided in VARARGIN. 
%
%   _Authors: Chung-Heng Yeh <chyeh@ee.columbia.edu>_
%
%   _Copyright 2012-2014 Chung-Heng Yeh_

function dx = mem_con_template(v, x, t, g, varargin)

    % Handle the optional input parameters. You can simply jump to line 59.
    % =====================================================================
    p = inputParser;
    % Add optinal argument using following format
    % addParamValue(p,'g',[0 0 0],@isnumeric);
    parse(p,varargin{:});
    
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
    dV = - g(1)*x(2).^3*x(3)*(v-50) - g(2)*x(1).^4*(v+77) - g(3)*(v+54.387);
    dx = [dV dn dm dh];    
end