function [DV DN] = hhn_vector_field(V,N,g)
    s = (max(max(V)) - min(min(V))) / (max(max(N)) - min(min(N)));
    H = 1-N;
    M =  0.1*(V+40) ./ (1-exp(-(V+40)/10)) / 4 ./ exp(-(V+65)/18);
    DN = zeros(size(V));
    a  = exp(-(V+55)/10)-1;
    DN( a==0 ) = ( 1-N( a==0 ) ) * 0.1 - N( a==0 ) .* (0.125.*exp(-(V( a==0 ) +65)/80));
    DN( a~=0 ) = (1-N( a~=0 )) .* (-0.01*(V(a~=0)+55)./a(a~=0)) - N( a~=0 ) .* (0.125*exp(-(V(a~=0)+65)/80));
    DN = DN*s;
    DV = (- g(1)*M.^3.*H.*(V-50) - g(2)*N.^4.*(V+77) - g(3)*(V+54.387)); 
    b = sqrt(DN.^2 + DV.^2);
    DN = DN./b/s;
    DV = DV./b;
end