function [out, J] = HH_Jacobian(x,g) 

    gNa = g(1);
    gK  = g(2);
    gL  = g(3);
    ENa = 50;
    EK  = -77;
    EL  = -54.387;
    v   = x(1);
    
    % Compute alpha and beta for n,m, and h.
    if v == -55
        a_n = 0.1;
    else
        a_n = 0.01*(v+55)/(1-exp(-(v+55)/10));
    end
    if v == -40
        a_m = 1;
    else
        a_m = 0.1*(v+40)/(1-exp(-(v+40)/10));
    end    
    a_h = (0.07*exp(-(v+65)/20));
    b_n =  0.125*exp(-(v+65)/80);
    b_m = 4*exp(-(v+65)/18);
    b_h = 1 / (exp(-(v+35)/10)+1);
    a_n_dv = 0.01*(1-(v+65)/10*exp(-(v+55)/10))/(1-exp(-(v+55)/10))^2;
    a_m_dv = 0.1*(1-(v+50)/10*exp(-(v+40)/10))/(1-exp(-(v+40)/10))^2;
    a_h_dv = -7/200*exp(-(v+65)/20);
    b_n_dv = -b_n/80;
    b_m_dv = -b_m/18;
    b_h_dv = 0.1*exp(-(v+35)/10)*b_h^2;
    
    % if only V is given, n, m, and h are set to the steady state.
    if numel(x) == 1,
        n = a_n/(a_n+b_n);
        m = a_m/(a_m+b_m);
        h = a_h/(a_h+b_h);
    else
        n = x(2);
        m = x(3);
        h = x(4);
    end
    
    
      
    
    % compute partial derivative of v'
    dDV_dV = (-gNa*m^3*h-gK*n^4-gL);
    dDV_dn = -4*gK*n^3*(v-EK);
    dDV_dm = -3*gNa*m^2*h*(v-ENa);
    dDV_dh = -gNa*m^3*(v-ENa);
   
    % compute partial derivative
    dDn_dn = - (a_n+b_n);
    dDm_dm = - (a_m+b_m);
    dDh_dh = - (a_h+b_h);
    
    dDn_dV = a_n_dv*(1-n) - b_n_dv*n;
    dDm_dV = a_m_dv*(1-m) - b_m_dv*m;
    dDh_dV = a_h_dv*(1-h) - b_h_dv*h;
    
    J = [dDV_dV dDV_dn dDV_dm dDV_dh ;...
         dDn_dV dDn_dn      0      0 ;...
         dDm_dV      0 dDm_dm      0 ;...
         dDh_dV      0      0 dDh_dh ;];
    
    out = eig(J);
    I = gNa*m^3*h*(v-ENa)+gK*n^4*(v-EK)+gL*(v-EL);
    %[I v n m h]
end
