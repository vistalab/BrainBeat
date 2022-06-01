%
% This script makes Supplemental Figure 6 from:
%
%%%%% Measuring brain beats: cardiac-aligned fast fMRI signals %%%%%
% Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell
%
%

clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%% Load respiratory data

all_phys = [];

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,1,1,1,1};

for ss = 1:length(sub_labels)    
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    rr = 1;    
    run_nr = run_nrs{ss}(rr);

    fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
    fmri_name = fullfile(dDir,fmri_BIDSname);
    ni = niftiRead(fmri_name);

    % get physio stuff we need:
    all_phys(ss).physio = physioCreate('nifti',ni);
    all_phys(ss).resp_onsets = physioGet(all_phys(ss).physio,'resp peaks');
end


%%
figure('Position',[0 0 800 800])

for ss = 1:length(sub_labels) 

    subplot(6,2,ss*2-1),hold on

    resp_srate = all_phys(ss).physio.resp.srate;
    resp_plot = zscore(all_phys(ss).physio.resp.data);
    plot((1:length(all_phys(ss).physio.resp.data))/resp_srate,resp_plot,'Color',[.2 .2 .2])
    plot(all_phys(ss).resp_onsets,resp_plot(round(all_phys(ss).resp_onsets*resp_srate)),'k*')
    xlim([0 all_phys(ss).resp_onsets(end)])
    xlabel('Time (s)')

    subplot(6,2,ss*2),hold on
    plot(all_phys(ss).resp_onsets(2:end),diff(all_phys(ss).resp_onsets),'Color',[.5 .5 .5])
    xlabel('Time (s)')
    xlim([0 all_phys(ss).resp_onsets(end)])

end

subplot(6,2,5),
ylabel('Respiratory signal (z-score)')

subplot(6,2,6),
ylabel('Breath-breath interval (s)')

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures','SuppFigureS6_respiration'))
print('-r300','-depsc2',fullfile(dDir,'derivatives','figures','SuppFigureS6_respiration'))


