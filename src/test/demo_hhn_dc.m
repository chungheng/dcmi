%%DEMO of HODGKIN_HUXLEY_dynamicClamp 
%   In this code, we demonstrate the usage of HODGKIN_HUXLEY_dynamicClamp. 
%   Each channel of Hodgkin-Huxley neuron is disabled by both 
%   pharmacological method as well as dynamic clamping
%
%   _Authors: Chung-Heng Yeh <chyeh@ee.columbia.edu>_
%
%   _Copyright 2012-2014 Chung-Heng Yeh_

% Initialize demo
addpath('../hhn/');
dt = 1e-5;
t = 0:dt:0.05;
I = ones(size(t))*40;

figure()
       
V = hodgkin_huxley_dynamicClamp(t,I);
subplot(4,1,1);
    plot(t,V(:,1));
    ylabel('Voltage, [mV]');title('HHN');xlim([t(1) t(end)]);

V = hodgkin_huxley_dynamicClamp(t,I,'Potassium','off');
U = hodgkin_huxley_dynamicClamp(t,I,'MemInitState',[0 0 1],...
                      'MEMCON',@(v,x,t) mem_con_template(v,x,t,[0 -36 0]));
subplot(4,1,2);
    plot(t,V(:,1),t,U(:,1));title('Potassium Channel Disabled');xlim([t(1) t(end)]);
    ylabel('Voltage, [mV]');legend('pharmacological','via Dynamic Clamp');
    
V = hodgkin_huxley_dynamicClamp(t,I,'Sodium','off');
U = hodgkin_huxley_dynamicClamp(t,I,'MemInitState',[0 0 1],...
                      'MEMCON',@(v,x,t) mem_con_template(v,x,t,[-120 0 0]));
subplot(4,1,3);
    plot(t,V(:,1),t,U(:,1));title('Sodium Channel Disabled');xlim([t(1) t(end)]);
    ylabel('Voltage, [mV]');legend('pharmacological','via Dynamic Clamp');

V = hodgkin_huxley_dynamicClamp(t,I,'Leaky','off');
U = hodgkin_huxley_dynamicClamp(t,I,'MemInitState',[0 0 1],...
                      'MEMCON',@(v,x,t) mem_con_template(v,x,t,[0 0 -0.3]));
subplot(4,1,4);
    plot(t,V(:,1),t,U(:,1));title('Leaky Channel Disabled');xlim([t(1) t(end)]);
    xlabel('Time, [s]');ylabel('Voltage, [mV]');
    legend('pharmacological','via Dynamic Clamp');


