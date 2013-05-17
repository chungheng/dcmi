%% 
clc;
dt = 1e-5;
t  = -0.1:dt:1;
I  = 0*ones(size(t));
gl = -1;

[V M] = hodgkin_huxley_BSDC(t,I,...
        'memcon', @mem_con,'MemInitState',[0 0 0],'MaxCon',[0 0 gl],'MaxConUpdate',[0 0 0.1]);
figure();
    subplot(2,1,1);
        plot(t*1e3,V(:,1));xlabel('time, [ms]');ylabel('voltage, [mV]');
    subplot(2,1,2);
        plot(t*1e3,M(:,3));xlabel('time, [ms]');ylabel('conductance, [mS]');