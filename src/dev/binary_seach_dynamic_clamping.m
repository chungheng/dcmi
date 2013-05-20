%% 
clc;
addpath('../hhn/');
%% Identification of leaky channel using binary search dynamic clamping
dt = 1e-5;
t  = -0.1:dt:2;
I  = 0*ones(size(t));

[V M S] = hodgkin_huxley_BSDC(t,I,'Verbose',true,...
          'memcon', @mem_con,'MemInitState',[0 0 0],'MaxCon',[0 0 -1.0],'MinCon',[0 0 0]);
%%
figure();
    subplot(3,1,1);
        plot(t*1e3,V(:,1));xlim(1e3*[t(1) t(end)]);ylabel('Membrane voltage, [mV]');
        ylim([-100 100]);
    subplot(3,1,2);
        plot(t*1e3,M(:,3));xlim(1e3*[t(1) t(end)]);ylim([-1.1 0.1]);
        hold on;plot( 1e3*[t(1) t(end)], [-0.3 -0.3],'--r');ylabel('Injected Conducatance');
        legend('Injected','Original','Location','SouthEast');
    subplot(3,1,3);
        plot(t*1e3,S(:,1),'g');xlim(1e3*[t(1) t(end)]);
        hold on;plot(t*1e3,S(:,2)-1.5,'k');xlim(1e3*[t(1) t(end)]);
        ylim([-2.1 1.6]);ylabel('Contral Signal');set(gca,'YTick',[]);xlabel('time, [ms]');
        legend('Reset Signal','Breakdown Indicator','Location','SouthEast');