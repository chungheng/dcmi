%% 
ss =  [-64.9964    0.3177    0.0530    0.5960];
g  =  [120 36 0.3];
[x fval flag] = fsolve( @(x) hhn(x,g), ss );

%% 
% Use Rinzel's model to reduce the model to 2 variable V and n, and plot the 
% corresponding vector field
g  =  [120 36 0.3];
v = -70:-60;
n = 0.25:.01:0.35;
[V, N] = meshgrid(v,n);

[DV DN] = hhn_vector_field(V, N, g);
figure();quiver(V,N,DV,DN,0.1);
%%
gL = -1:.05:1;
X  = zeros(numel(gL), 4);
for i = 1:numel(gL)
    [x fval flag] = fsolve( @(x) hhn(x,[120 36 gL(i)]), ss );
    if flag == 1
        X(i,:)  = x;
    end
end

figure();plot3(X(:,1),X(:,2),gL);grid on;
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

