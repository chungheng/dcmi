%% 
ss =  [-64.9964    0.3177    0.0530    0.5960];
g  =  [120 36 0.3];

g_hat = -1:0.1:1;
E = zeros(numel(g_hat),4);

for i = 1:numel(g_hat)
    E(i,:) = HH_Jacobian(ss,[120 36 g_hat(i)]);
end
figure();cs = colormap(hsv(7));
for i = 1:4
    hold on;plot(real(E(:,i)),imag(E(:,i)),'color',cs(i,:));
end

