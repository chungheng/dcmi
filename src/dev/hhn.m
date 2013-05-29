function dx = hhn( x, g )
    %g = [120 36 0.3];
    v = x(1);
    a = exp(-(v+55)/10)-1;
    if a == 0
        dn = (1-x(2)) * 0.1 - x(2) * (0.125*exp(-(v+65)/80));
    else
        dn = (1-x(2)) * (-0.01*(v+55)/a) - x(2) * (0.125*exp(-(v+65)/80));
    end

    a = exp(-(v+40)/10)-1;
    if a == 0
        dm = (1-x(3)) - x(3) * (4*exp(-(v+65)/18));
    else
        dm = (1-x(3)) * (-0.1*(v+40)/a) - x(3) * (4*exp(-(v+65)/18));
    end
    dh = (1-x(4)) * (0.07*exp(-(v+65)/20)) - x(4) / (exp(-(v+35)/10)+1);
    I_nmh = - g(1)*x(3).^3*x(4)*(v-50) - g(2)*x(2).^4*(v+77) - g(3)*(v+54.387);
    dx    = [I_nmh dn dm dh]; 
end

