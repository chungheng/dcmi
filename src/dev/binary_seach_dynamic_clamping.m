%% 
clc;
addpath('../hhn/');
%% Identification of leaky channel using binary search dynamic clamping
dt = 1e-5;
t  = -0.1:dt:2;
I  = 0*ones(size(t));
gi = -1;
g  = -0.3;
[V M S] = hodgkin_huxley_BSDC(t,I,'Verbose',true,...
          'memcon', @mem_con,'MemInitState',[0 0 0],'MaxCon',[0 0 gi],'MinCon',[0 0 0]);
%%
fig_lc = figure('Name','Limit cycle');
    set(gca,'nextplot','replacechildren');
    plot(V(:,1),V(:,2));xlim([-120 50]);ylim([0 1]);
    hold on; anim1 = plot(V(1,1),V(1,2),'or');set(anim1,'EraseMode','normal');

fig_cp = figure('Name','Control Panel');
    subplot(3,1,1);set(gca,'nextplot','replacechildren');
        plot(t*1e3,V(:,1));xlim(1e3*[t(1) t(end)]);ylabel('Membrane voltage, [mV]');
        ylim([-100 100]);
    subplot(3,1,2);set(gca,'nextplot','replacechildren');
        plot(t*1e3,M(:,3));xlim(1e3*[t(1) t(end)]);ylim([gi-0.1 0.1]);
        hold on;plot( 1e3*[t(1) t(end)], [g g],'--r');ylabel('Injected Conducatance');
        legend('Injected','Original','Location','SouthEast');
    subplot(3,1,3);set(gca,'nextplot','replacechildren');
        plot(t*1e3,S(:,1),'g');xlim(1e3*[t(1) t(end)]);
        hold on;plot(t*1e3,S(:,2)-1.5,'k');xlim(1e3*[t(1) t(end)]);
        ylim([-2.1 1.6]);ylabel('Contral Signal');set(gca,'YTick',[]);xlabel('time, [ms]');
        legend('Reset Signal','Breakdown Indicator','Location','SouthEast');
    subplot(3,1,1);hold on;anim2 = plot(t(1)*1e3,V(1,1),'ro');set(anim2,'EraseMode','normal');
    subplot(3,1,2);hold on;anim3 = plot(t(1)*1e3,M(1,3),'ro');set(anim3,'EraseMode','normal');
    subplot(3,1,3);hold on;anim4 = plot(t(1)*1e3,S(1,1),'ro');set(anim4,'EraseMode','normal');
%%
for i = 1:75:numel(t)
    set(anim1,'XData',V(i,1),'YData',V(i,2));
    set(anim2,'XData',t(i)*1e3,'YData',V(i,1));
    set(anim3,'XData',t(i)*1e3,'YData',M(i,3));
    set(anim4,'XData',t(i)*1e3,'YData',S(i,1));
    pause(dt/5);
end
        
        %%
%{
vidObj = VideoWriter('bonus.avi');
vidObj.FrameRate = 30;
vidObj.Quality = 100;
open(vidObj);
%}
%%

%% Identification of potassium channel using binary search dynamic clamping
dt = 1e-5;
t  = -0.1:dt:2;
I  = 0*ones(size(t));
gi = -60;
g  = -36;
[V M S] = hodgkin_huxley_BSDC(t,I,'Verbose',true,...
          'memcon', @mem_con,'MemInitState',[0 0 0],'MaxCon',[0 gi 0],'MinCon',[0 0 0]);
figure();
    subplot(3,1,1);
        plot(t*1e3,V(:,1));xlim(1e3*[t(1) t(end)]);ylabel('Membrane voltage, [mV]');
        ylim([-100 100]);
    subplot(3,1,2);
        plot(t*1e3,M(:,2));xlim(1e3*[t(1) t(end)]);ylim([gi-0.1 0.1]);
        hold on;plot( 1e3*[t(1) t(end)], [g g],'--r');ylabel('Injected Conducatance');
        legend('Injected','Original','Location','SouthEast');
    subplot(3,1,3);
        plot(t*1e3,S(:,1),'g');xlim(1e3*[t(1) t(end)]);
        hold on;plot(t*1e3,S(:,2)-1.5,'k');xlim(1e3*[t(1) t(end)]);
        ylim([-2.1 1.6]);ylabel('Contral Signal');set(gca,'YTick',[]);xlabel('time, [ms]');
        legend('Reset Signal','Breakdown Indicator','Location','SouthEast');
        
%% Identification of sodium channel using binary search dynamic clamping
dt = 1e-5;
t  = -0.1:dt:3;
I  = 0*ones(size(t));
gi = -140;
g  = -120;
[V M S] = hodgkin_huxley_BSDC_sodium(t,I,'Verbose',true,...
          'memcon', @mem_con,'MemInitState',[0 0 0],'MaxCon',[gi -36 0],'MinCon',[0 -36 0]);
figure();
    subplot(3,1,1);
        plot(t*1e3,V(:,1));xlim(1e3*[t(1) t(end)]);ylabel('Membrane voltage, [mV]');
        ylim([-100 100]);
    subplot(3,1,2);
        plot(t*1e3,M(:,1));xlim(1e3*[t(1) t(end)]);ylim([gi-0.1 0.1]);
        hold on;plot( 1e3*[t(1) t(end)], [g g],'--r');ylabel('Injected Conducatance');
        legend('Injected','Original','Location','SouthEast');
    subplot(3,1,3);
        plot(t*1e3,S(:,1),'g');xlim(1e3*[t(1) t(end)]);
        hold on;plot(t*1e3,S(:,2)-1.5,'k');xlim(1e3*[t(1) t(end)]);
        ylim([-2.1 1.6]);ylabel('Contral Signal');set(gca,'YTick',[]);xlabel('time, [ms]');
        legend('Reset Signal','Decreasing Indicator','Location','SouthEast');