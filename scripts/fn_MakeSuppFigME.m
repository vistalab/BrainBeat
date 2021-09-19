function [] = fn_MakeSuppFigME(d1,d2,lns0,rr,t2s,srate,ppg_onsets)

%% Make sure we plot TE/T2s as well, this is the slope we are calculating.

x_l = [90 100];
figure('Position',[0 0 600 800])
tt = (1:size(d1,3))/srate;
subplot(5,2,1), hold on
plot(tt,squeeze(d1(27,50,:)),'b')
plot(ppg_onsets,mean(d1(27,50,:)),'k.')    
ylabel('Ln(S_T_E_1)'),xlabel('time(s)')
xlim(x_l)

subplot(5,2,3), hold on
plot(tt,squeeze(d2(27,50,:)),'g')
plot(ppg_onsets,mean(d2(27,50,:)),'k.')    
ylabel('Ln(S_T_E_1)'),xlabel('time(s)')
xlim(x_l)

subplot(5,2,[2 4]), hold on
plot(tt,squeeze(d1(27,50,:)),'b')
plot(tt,squeeze(d2(27,50,:)),'g')
plot([tt; tt],[squeeze(d1(27,50,:)) squeeze(d2(27,50,:))]')
plot(ppg_onsets,mean(d1(27,50,:)),'k.')    
plot(ppg_onsets,mean(d2(27,50,:)),'k.')    
% plot(ppg_onsets,.5*[mean(d2(27,50,:))+mean(d1(27,50,:))],'k.')    
xlim(x_l),ylim([8 9])

subplot(5,2,6), hold on
plot(tt,squeeze(lns0(27,50,:)),'r')
plot(ppg_onsets,mean(lns0(27,50,:)),'k.')
ylabel('Ln(S_0)')
xlim(x_l)

subplot(5,2,8), hold on, ylabel('rr')
plot(tt,squeeze(rr(27,50,:)),'m')
plot(ppg_onsets,mean(rr(27,50,:)),'k.')
xlim(x_l)

subplot(5,2,10), hold on, ylabel('T2s')
plot(tt,squeeze(t2s(27,50,:)),'m')
plot(ppg_onsets,mean(t2s(27,50,:)),'k.')
xlim(x_l)
xlabel('time(s)')

figName = fullfile(bbPath,'local','ME_ExampleFig');

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',figName)
print('-painters','-r300','-depsc',figName)

