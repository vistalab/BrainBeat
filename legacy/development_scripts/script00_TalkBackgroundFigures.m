clear all
close all
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

figure('Position',[0 0 400 200])

RT = .5; 
hrf = spm_hrf(RT); 
plot(0:RT:32, hrf,'k','LineWidth',2);
xlim([0 32]),ylim([min(hrf) max(hrf)])
xlabel('time (s)')
ylabel('')
set(gca,'XTick',[0:10:30],'YTick',[],'FontSize',20)
box off

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/methods/HRF'])
