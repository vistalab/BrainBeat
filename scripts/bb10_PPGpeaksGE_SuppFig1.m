%
% This script makes Supplemental Figure 1 from:
%
%%%%% Measuring brain beats: cardiac-aligned fast fMRI signals %%%%%
% Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell
%
%

clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%%
sub_label = 'sub-5';
ses_label = 'ses-1';

% load fMRI
ni = niftiRead(fullfile(dDir,sub_label,ses_label,'func',...
    [sub_label '_' ses_label '_task-rest_acq-4mmFA48_run-1_bold.nii.gz']));
% get our physio analysis
physio = physioCreate('nifti',ni);
ppg_onsets = physioGet(physio,'ppg peaks');

% get scanner data and triggers
ppgData = bids.util.tsvread(fullfile(dDir,'sourcedata','sub-5',...
    [sub_label '_ses-1_task-rest_acq-4mmFA48_run-1_PPGData_raw.tsv']));
ppgTrig = bids.util.tsvread(fullfile(dDir,'sourcedata','sub-5',...
    [sub_label '_ses-1_task-rest_acq-4mmFA48_run-1_PPGTrig_raw.tsv']));


%%

figure('Position',[0 0 600 200]),hold on
% subplot(2,1,1),hold on
ppg_onsets_samples = floor(ppg_onsets*physio.ppg.srate);

% offset due to start scan after PPG start
ppg_offset = (length(physio.ppg.rawdata)-length(physio.ppg.data))/physio.ppg.srate;

tt = ppg_offset + (1:length(physio.ppg.data))/physio.ppg.srate;

plot(tt,physio.ppg.data,'k')
plot(ppg_onsets+ppg_offset,physio.ppg.data(ppg_onsets_samples),'g*')
xlim([230 250])

plot(ppgTrig/physio.ppg.srate,ppgData(ppgTrig),'r*')
xlim([230 250])

xlabel('Time (s)')
ylabel('PPG')

% set(gcf,'PaperPositionMode','auto')
print('-depsc2','-r300',fullfile(dDir,'derivatives','figures','SuppFigureS1_PPGpeaksGE'))
print('-dpng','-r300','-painters',fullfile(dDir,'derivatives','figures','SuppFigureS1_PPGpeaksGE'))
    
    

