%% Initialization
clc; close all; clear all;
% Add path 
addpath('../hhn/');

%% Counter-balance the leaky channel
% The maximum channel conductance of the injected memconductance is initialize 
% to a zero vector in the beginning of simulation, and set to the MAXCON 
% after RESETTIME. 
dt = 1e-5;
t  = -0.1:dt:.3;
I  = 0*ones(size(t));
gl = 0.1:0.1:1;
le = cell(1,numel(gl));
figure(); cs = colormap(hsv(numel(gl)+4));
for i = 1:numel(gl)
    V  = hodgkin_huxley_dynamicClamp(t,I,'ResetTime',0,...
         'memcon', @mem_con,'MemInitState',[0 0 0],'MaxCon',[0 0 -gl(i)]);
    hold on;plot(t*1e3,V(:,1),'Color',cs(i,:));
    le{i} = num2str(-gl(i));
end
ylim([-100 20]);legend(le);xlim([-50 100]);
xlabel('time, [ms]');ylabel('voltage, [mV]');

%% Counter-balance the Sodium Channel channel
% The maximum channel conductance of the injected memconductance is initialize 
% to a zero vector in the beginning of simulation, and set to the MAXCON 
% after RESETTIME. 
dt  = 1e-5;
t   = -0.1:dt:.3;
I   = 5*ones(size(t));
gNa = 115:125;
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
gK = 30:40;
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