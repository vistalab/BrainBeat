clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';
% chdir(dDir)

%% Plot fMRI signals from DKT segmentation in one subject

sub_labels = {'1'}; 
ses_labels = {'1'}; 
acq_labels = {'4mmFA48'};
run_nrs = {[1]};

ss = 1;%:length(sub_labels) % subjects/ses/acq
rr = 1;% run_nr
sub_label = sub_labels{ss};
ses_label = ses_labels{ss};
acq_label = acq_labels{ss};

run_nr = run_nrs{ss}(rr);
        
dkt_table = readtable('dkt_areas.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dktNames = dkt_table.label;
dktCodes = dkt_table.label_nr;

fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
fmri_name = fullfile(dDir,fmri_BIDSname);
ni = niftiRead(fmri_name);

% get physio stuff we need:
physio      = physioCreate('nifti',ni);
ppg_onsets  = physioGet(physio,'ppg peaks');
ppgCurve    = physioGet(physio,'ppg ppgcurve');
ppgCurveT   = physioGet(physio,'ppg ppgtcurve');
srate       = 1/bbGet(ni,'tr');

save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

ppgResp = niftiRead([save_name_base '_PPGtrigResponse.nii.gz']); % ppg triggered time series
ppgT = load([save_name_base '_PPGtrigResponseT.mat'],'t');

%%
niSegm = niftiRead([save_name_base '_combineSegm.nii.gz']);
roiNames = {'GM','WM','Ventricles','CSF','Veno'};
% Freesurfer for GM, WM, Ventricles, CSF from SPM and Venogram

% Voxel vector with segmentation labels 1:5 for [GM WM Ventricles CSF Veno]
segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

% load DKT atlas
niDKT = niftiRead([save_name_base '_r_DKTatlas_aseg.nii.gz']);
% matrix to vector
segmVectDKT = reshape(niDKT.data,[size(niDKT.data,1) * size(niDKT.data,2) * size(niDKT.data,3)],1);

%%%% put the PPG triggered data in a matrix Voxel X Time
respMat = reshape(ppgResp.data,[size(ppgResp.data,1) * size(ppgResp.data,2) * size(ppgResp.data,3)],size(ppgResp.data,4));

% initialize average PPG triggered response across ROIs
avResp = zeros(size(ppgResp.data,4),length(roiNames)+1);
for kk = 1:length(roiNames)
    avResp(:,kk) = mean(respMat(segmVect==kk,:),1);
end

this_areas = 1; % cingulate
avResp(:,kk+1) = mean(respMat(ismember(segmVectDKT,dktCodes(this_areas)),:),1);


%%%% put the fMRI data in a matrix Voxel X Time
sigMat = reshape(ni.data,[size(ni.data,1) * size(ni.data,2) * size(ni.data,3)],size(ni.data,4));

% initialize average total time signal
avSig = zeros(size(ni.data,4),length(roiNames)+1);
for kk = 1:length(roiNames)
    avSig(:,kk) = mean(sigMat(segmVect==kk,:),1);
end
avSig(:,kk+1) = mean(sigMat(ismember(segmVectDKT,dktCodes(this_areas)),:),1);
avSig(1:5,:) = NaN; % first scans to NaN

roiNames(length(roiNames)+1) = dktNames(this_areas);

%%

areas_plot = [6 3 5];

tt = [1:size(avSig,1)]/srate;

figure('Position',[0 0 500 400])
for kk = 1:length(areas_plot)
    thisArea = areas_plot(kk);
    
    sig_plot = avSig(:,thisArea)';
    % detrend
    x = 6:length(sig_plot);
    p = polyfit(x,sig_plot(x),1);    

    % percent modulation
    mean_factor = mean(sig_plot(x)); % mean
    sig_plot = 100*(sig_plot - (p(1)*[1:length(sig_plot)] + p(2)))./mean_factor;
    
    subplot(length(areas_plot)+1,2,kk*2-1),hold on
    plot(tt,sig_plot,'.','Color',[0 .6 .8],'MarkerSize',10)
    plot(tt,sig_plot,'Color',[0 .6 .8],'LineWidth',1)
%     plot([ppg_onsets ppg_onsets],[nanmean(avSig(:,thisArea))*100 nanmean(avSig(:,thisArea))*100],'k')
    plot([ppg_onsets ppg_onsets],[min(avResp(:,thisArea)*100) max(avResp(:,thisArea)*100)],'k')
    xlim([85 95]) 
end

for kk = 1:length(areas_plot)
    thisArea = areas_plot(kk);
    subplot(length(areas_plot)+1,2,kk*2),hold on
    plot(ppgT.t,avResp(:,thisArea)*100,'Color',[0 .6 .8],'LineWidth',2)
    plot([0 0],[min(avResp(:,thisArea)*100) max(avResp(:,thisArea)*100)],'k')
    xlim([min(ppgT.t) max(ppgT.t)])
    title(roiNames{thisArea})
end

% ppg for completeness
subplot(length(areas_plot)+1,2,length(areas_plot)*2+1),hold on
plot([ppg_onsets ppg_onsets],[-50 100],'k')
plot([1:length(physio.ppg.data)]/physio.ppg.srate,physio.ppg.data)
xlim([85 95])

subplot(length(areas_plot)+1,2,length(areas_plot)*2+2),hold on
plot(ppgCurveT,ppgCurve,'Color',[0 .6 .8],'LineWidth',2)
plot([0 0],[min(ppgCurve) max(ppgCurve)],'k')
xlim([min(ppgCurveT) max(ppgCurveT)])
title('ppg')

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir '/figures/segmentation/s' int2str(s_nr) '_scan' int2str(scan_nr) '_tracessegm'])
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/s' int2str(s_nr) '_scan' int2str(scan_nr) '_tracessegm'])


