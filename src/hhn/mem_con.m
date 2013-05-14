% Define Default Hodgkin-Huxley Nueron ODEs. This function handle is 
% written by Yiyin Zhou.
function dx =mem_con(v,x,g)
    %g = [0 0 -0.4];
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