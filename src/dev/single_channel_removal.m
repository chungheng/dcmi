%% Initialization
clc; close all; clear all;
% Add path 
addpath('../hhn/');
%% Reset process
% Demonstrate the dynamics of the state varibles of both mem-conductance 
% and the neuron model during RESETTIME.
dt = 1e-5;
t  = -0.05:dt:.3;
I  = 0*ones(size(t));

[V M] = hodgkin_huxley_BSDC(t,I,'ResetPeriod',0,...
         'memcon', @mem_con,'MemInitState',[0 0 0],'MaxCon',[0 0 0],'MinCon',[0 0 0]);

figure();plot(t*1e3,V(:,1));xlim([-50 0]);ylim([-80 20]);
         xlabel('time, [ms]');ylabel('voltage, [mV]');
figure();plot(t*1e3,V(:,2),'r',t*1e3,M(:,4),'--r',...
              t*1e3,V(:,3),'g',t*1e3,M(:,5),'--g',...
              t*1e3,V(:,4),'b',t*1e3,M(:,6),'--b');xlim([-50 0]);
         xlabel('time, [ms]');ylim([0 1]);ylabel('State Variable');
         legend('neu-n','mem-n','neu-m','mem-m','neu-h','mem-h');
%% Counter-balance the leaky channel
% The maximum channel conductance of the injected memconductance is initialize 
% to a zero vector in the beginning of simulation, and set to the MAXCON 
% after RESETTIME. 
dt = 1e-5;
t  = -0.1:dt:.3;
I  = 0*ones(size(t));
gl = 0.25:0.01:0.35;
le = cell(1,numel(gl));
figure(); cs = colormap(hsv(numel(gl)+4));
for i = 1:numel(gl)
    V  = hodgkin_huxley_dynamicClamp(t,I,'ResetTime',0,...
         'memcon', @mem_con,'MemInitState',[0 0 0],'MaxCon',[0 0 -gl(i)]);
    hold on;plot(t*1e3,V(:,1),'Color',cs(i,:));
    le{i} = num2str(-gl(i));
end
ylim([-150 20]);legend(le);xlim([-50 100]);
xlabel('time, [ms]');ylabel('voltage, [mV]');

%% Counter-balance the Sodium Channel channel
% The maximum channel conductance of the injected memconductance is initialize 
% to a zero vector in the beginning of simulation, and set to the MAXCON 
% after RESETTIME. 
dt  = 1e-5;
t   = -0.1:dt:.3;
I   = 5*ones(size(t));
gNa = 119.5:0.1:120.5;
le  = cell(1,numel(gNa));
figure(); cs = colormap(hsv(numel(gNa)+4));
for i = 1:numel(gNa)
    V  = hodgkin_huxley_dynamicClamp(t,I,'ResetTime',0,...
         'memcon', @mem_con,'MemInitState',[0 0 0],'MaxCon',[-gNa(i) -36 0]);
    hold on;plot(t*1e3,V(:,1),'Color',cs(i,:));
    le{i} = num2str(-gNa(i));
end
ylim([-80 0]);legend(le);xlim([-50 100]);
xlabel('time, [ms]');ylabel('voltage, [mV]');

%% Counter-balance the Potassium Channel channel
% The maximum channel conductance of the injected memconductance is initialize 
% to a zero vector in the beginning of simulation, and set to the MAXCON 
% after RESETTIME. 
dt = 1e-5;
t  = -0.1:dt:.3;
I  = 0*ones(size(t));
gK = 35.5:0.1:36.5;
le = cell(1,numel(gK));
figure(); cs = colormap(hsv(numel(gK)+4));
for i = 1:numel(gK)
    V  = hodgkin_huxley_dynamicClamp(t,I,'ResetTime',0,...
         'memcon', @mem_con,'MemInitState',[0 0 0],'MaxCon',[0 -gK(i) 0]);
    hold on;plot(t*1e3,V(:,1),'Color',cs(i,:));
    le{i} = num2str(-gK(i));
end
ylim([-100 100]);legend(le);xlim([-50 100]);
xlabel('time, [ms]');ylabel('voltage, [mV]');