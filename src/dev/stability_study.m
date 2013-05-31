%% Stability study of pseudo leaky channel injection
clc;
addpath('../hhn/');
%%
% Plot the null-solution for gL of different value. There is no solution 
% when gL is negative
ss =  [-64.9964    0.3177    0.0530    0.5960];
gL = -1:.05:1;
X  = zeros(numel(gL), 4);
f  = zeros(1,numel(gL));
for i = 1:numel(gL)
    [x fval flag] = fsolve( @(x) hhn(x,[120 36 gL(i)]), ss );
    X(i,:)  = x;
    f(i) = flag;
end
idx = find(f == 1);

le = cell(1,numel(idx)-1);
figure();cs = colormap( hsv( numel(idx)+1 ) );

for i = 1:numel(idx)
    hold on;plot(X(idx(i),1),X(idx(i),2),'o','color',cs(i,:));
    le{i} = num2str(gL(idx(i)));
end
grid on;axis([-100 -55 0 0.5]);legend(le,'Location','NorthWest');
xlabel('Voltage, [mV]');ylabel('N');
%%
% Use Rinzel's model to reduce the model to 2 variable V and n, and plot 
% the corresponding vector field. The vector field for different value of
% gL is cascaded into a gif movie.
v = -100:2:-60;
n = 0:0.02:0.5;
[V, N] = meshgrid(v,n);

gL = -1:.05:1;

figure();
filename = 'gl_vf.gif';
for i = 1:numel(gL)
    set(gcf,'Color','W');
    [DV DN] = hhn_vector_field(V, N, [120 36 gL(i)]);
    plot(X(idx,1),X(idx,2),'k');
    hold on;plot( ss(1), ss(2), 'go','Markerfacecolor','g');
    if any(idx == i)
        hold on;plot( X(i,1), X(i,2), 'ro','Markerfacecolor','r');
        le = {'Trace of Steady State',...
              'Steady state of the oringinal system',...
              'Steady state of the new system' };
    else
        le = {'Trace of Steady State','Steady state of the oringinal system'};
    end
    
    quiver(V,N,DV,DN,0.5,'Marker','o','ShowArrowHead','Off');
    xlim([-105 -55]);ylim([0 0.5]);title(['g*_L=',num2str(gL(i),'%3.2f')]);
    legend(le);
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
        clf
    end    
end
%%
% Print and plot trajectory of the eigenvalue of the Jacobian matrix
ss = [-64.9964    0.3177    0.0530    0.5960];
gL = [0 0.0022127 .1:.1:1];
E  = zeros(4,numel(gL));
for i = 1:numel(gL)
    x = fsolve( @(x) hhn(x,[120 36 gL(i)]), ss );
    E(:,i) = HH_Jacobian(x,[120 36 gL(i)]);
end
figure();
for i = 1:numel(gL)
    fprintf( '%1.4f & %1.4f & %1.4f%+1.4fi & %1.4f%+1.4fi & %1.4f\n', ...
            gL(i), E(1,i), real(E(2,i)), imag(E(2,i)), real(E(3,i)), imag(E(3,i)), E(4,i) );
end
plot( [-8 2], [0 0], '-k', [0 0], [-0.8 0.8], 'k');
hold on;plot( real(E(1,:)),imag(E(1,:)),'r','linewidth',1.5);
hold on;plot( real(E(2,:)),imag(E(2,:)),'g','linewidth',1.5);
hold on;plot( real(E(3,:)),imag(E(3,:)),'m','linewidth',1.5);
hold on;plot( real(E(4,:)),imag(E(4,:)),'b','linewidth',1.5);
xlim([-8 2]);ylim([-0.8 0.8]);
%% 
% Create a animation of eigenvalue trajectory
ss = [-64.9964    0.3177    0.0530    0.5960];
gL = [0:0.0005:0.00221 0.002:0.008:.1 .1:.05:1];
E  = zeros(4,numel(gL));
for i = 1:numel(gL)
    x = fsolve( @(x) hhn(x,[120 36 gL(i)]), ss );
    E(:,i) = HH_Jacobian(x,[120 36 gL(i)]);
end
figure('Color','w');clf
filename = 'gl_eig.gif';
plot( [-8 2], [0 0], '-k', [0 0], [-0.8 0.8], 'k');
hold on; plot( real(E(1,:)),imag(E(1,:)),'r');
hold on; plot( real(E(2,:)),imag(E(2,:)),'g');
hold on; plot( real(E(3,:)),imag(E(3,:)),'m');
hold on; plot( real(E(4,:)),imag(E(4,:)),'b');
hold on; h1 = plot(real(E(1,i)),imag(E(1,i)),'or','MarkerFaceColor','r');
hold on; h2 = plot(real(E(2,i)),imag(E(2,i)),'og','MarkerFaceColor','g');
hold on; h3 = plot(real(E(3,i)),imag(E(3,i)),'om','MarkerFaceColor','m');
hold on; h4 = plot(real(E(4,i)),imag(E(4,i)),'ob','MarkerFaceColor','b');
for i = 1:numel(gL)
    set(h1,'XData',real(E(1,i)),'YData',imag(E(1,i)));
    set(h2,'XData',real(E(2,i)),'YData',imag(E(2,i)));
    set(h3,'XData',real(E(3,i)),'YData',imag(E(3,i)));
    set(h4,'XData',real(E(4,i)),'YData',imag(E(4,i)));
    title(['g_L=',num2str(gL(i),'%3.2f')]);
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
    end  
end
xlim([-8 2]);ylim([-0.8 0.8]);
%% 
% Highlight the vector field 
v = -100:2:-60;
n = 0:0.02:0.5;
[V, N] = meshgrid(v,n);
% gL = 0.6
gl = 0.6;
[DV DN] = hhn_vector_field(V, N, [120 36 gl]);
idx = find( abs(gL-gl)<0.01, 1, 'first');
figure();plot( X(idx,1) , X(idx,2) ,'ro',ss(1),ss(2),'go');
         hold on;quiver(V,N,DV,DN,0.5,'Marker','o','ShowArrowHead','Off');
         xlim([-105 -55]);ylim([0 0.5]);title(['g*_L=',num2str(gl,'%3.2f')]);
         legend('Steady state of the new system','Steady state of the oringinal system');
% gL = 0.3
gl = 0.3;
[DV DN] = hhn_vector_field(V, N, [120 36 gl]);
idx = find( abs(gL-gl)<0.01, 1, 'first');
figure();plot( X(idx,1) , X(idx,2) ,'ro',ss(1),ss(2),'go');
         hold on;quiver(V,N,DV,DN,0.5,'Marker','o','ShowArrowHead','Off');
         xlim([-105 -55]);ylim([0 0.5]);title(['g*_L=',num2str(gl,'%3.2f')]);
         legend('Steady state of the new system','Steady state of the oringinal system');
% gL = 0
gl = 0.0;
[DV DN] = hhn_vector_field(V, N, [120 36 gl]);
idx = find( abs(gL-gl)<0.01, 1, 'first');
figure();plot( X(idx,1) , X(idx,2) ,'ro',ss(1),ss(2),'go');
         hold on;quiver(V,N,DV,DN,0.5,'Marker','o','ShowArrowHead','Off');
         xlim([-105 -55]);ylim([0 0.5]);title(['g*_L=',num2str(gl,'%3.2f')]);
         legend('Steady state of the new system','Steady state of the oringinal system');
% gL = 0.3
gl = -0.3;
[DV DN] = hhn_vector_field(V, N, [120 36 gl]);
figure();plot(ss(1),ss(2),'go');hold on;
         quiver(V,N,DV,DN,0.5,'Marker','o','ShowArrowHead','Off');
         %hold on;plot( X(idx,1) , X(idx,2) , 'ro');
         xlim([-105 -55]);ylim([0 0.5]);title(['g*_L=',num2str(gl,'%3.2f')]);
         legend('Steady state of the oringinal system');

%% Stability study of pseudo Potassium channel injection
%%
% Plot the null-solution for gK of different value. There is no solution 
% when gL is negative
ss =  [40    0.3177    0.0530    0.5960];
gK = -5:5;
X  = zeros(numel(gK), 4);
f  = zeros(1,numel(gK));
for i = 1:numel(gK)
    options = optimset('MaxFunEvals', 5000);
    [x fval flag] = fsolve( @(x) hhn(x,[120 gK(i) 0.3]), ss, options  );
    if flag == 1
        HH_Jacobian(x,[120 gK(i) 0.3])
    end
    X(i,:)  = x;
    f(i) = flag;
end
idx = find(f == 1);

le = cell(1,numel(idx)-1);
figure();cs = colormap( hsv( numel(idx)+1 ) );

for i = 1:numel(idx)
    hold on;plot(X(idx(i),1),X(idx(i),2),'o','color',cs(i,:));
    le{i} = num2str(gK(idx(i)));
end
grid on;axis([-80 80 0 1]);legend(le,'Location','NorthWest');
xlabel('Voltage, [mV]');ylabel('N');

%% 
% Highlight the vector field 
ss =  [-64.9964    0.3177    0.0530    0.5960];
v = -80:10:80;
n = 0:0.02:1;
[V, N] = meshgrid(v,n);
% gL = 0.6
gk = 2;
[DV DN] = hhn_vector_field(V, N, [120 gk 0.3]);
idx = find( abs(gK-gk)<0.01, 1, 'first');
figure();plot( X(idx,1) , X(idx,2) ,'ro',ss(1),ss(2),'go');
         hold on;quiver(V,N,DV,DN,0.5,'Marker','o','ShowArrowHead','Off');
         xlim([-80 80]);ylim([0 1]);title(['g*_L=',num2str(gk,'%3.2f')]);
         legend('Steady state of the new system','Steady state of the oringinal system');
% gL = 0.3
gk = 0;
[DV DN] = hhn_vector_field(V, N, [120 gk 0.3]);
idx = find( abs(gK-gk)<0.01, 1, 'first');
figure();plot( X(idx,1) , X(idx,2) ,'ro',ss(1),ss(2),'go');
         hold on;quiver(V,N,DV,DN,0.5,'Marker','o','ShowArrowHead','Off');
         xlim([-80 80]);ylim([0 1]);title(['g*_L=',num2str(gk,'%3.2f')]);
         legend('Steady state of the new system','Steady state of the oringinal system');
% gL = 0
gk = -2;
[DV DN] = hhn_vector_field(V, N, [120 gk 0.3]);
idx = find( abs(gK-gk)<0.01, 1, 'first');
figure();plot( X(idx,1) , X(idx,2) ,'ro',ss(1),ss(2),'go');
         hold on;quiver(V,N,DV,DN,0.5,'Marker','o','ShowArrowHead','Off');
         xlim([-80 80]);ylim([0 1]);title(['g*_L=',num2str(gk,'%3.2f')]);
         legend('Steady state of the new system','Steady state of the oringinal system');
% gL = 0.3
gk = -4;
[DV DN] = hhn_vector_field(V, N, [120 gk 0.3]);
idx = find( abs(gK-gk)<0.01, 1, 'first');
figure();plot( X(idx,1) , X(idx,2) ,'ro',ss(1),ss(2),'go');
         hold on;quiver(V,N,DV,DN,0.5,'Marker','o','ShowArrowHead','Off');
         xlim([-80 80]);ylim([0 1]);title(['g*_L=',num2str(gk,'%3.2f')]);
         legend('Steady state of the new system','Steady state of the oringinal system');