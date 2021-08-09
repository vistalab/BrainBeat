clear all
close all

% script produced transformation matrix between the anatomy and functional
% scans, walk through step-by-step

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% chdir(dDir)

%%
%% Overlay anat and reliability (COD)
%%

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
        
% Get the anatomicals:
t1w_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'anat',...
            ['sub-' sub_label '_ses-' ses_label '_T1w.nii.gz']);
niAnatomy = niftiRead(fullfile(dDir,t1w_BIDSname));

save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

% % Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));

% load coregistration matrix (for the functionals):
load([save_name_base '_AcpcXform_new.mat']);
acpcXform = acpcXform_new; clear acpcXform_new

%%%% Overlay 1: cod and anatomy
ppgRname = [save_name_base '_codPPG.nii.gz'];
ppgR = niftiRead(ppgRname); % correlation with PPG

%% Subject 1 Sagittal slices T1 + COD
curPos = [-10,1,-20]; 
sliceThisDim = 1; 
imDims=[-90 -120 -120; 90 130 90];
overlayPlot = ppgR;
cod_th = 0.8;

for kk = -60:10:60
    curPos(1) = kk;
    bbOverlayFuncAnat(overlayPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,1,1,cod_th);
    title(['slice ' int2str(kk) ', R^2>' num2str(cod_th,3)])
    set(gcf,'PaperPositionMode','auto')    
%     print('-painters','-r300','-dpng',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_orient' int2str(sliceThisDim) '_slice' int2str(kk)])
end

%% Subject 1 Axial slices T1 + COD
curPos = [-10,1,-20]; 
sliceThisDim = 3; 
imDims=[-90 -120 -120; 90 130 90];
overlayPlot = ppgR;
cod_th = 0.8;
for kk = -85:5:25
    curPos(3) = kk;
    bbOverlayFuncAnat(overlayPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,1,1,cod_th);
    title(['slice ' int2str(kk) ', R^2>' num2str(cod_th,3)])
    set(gcf,'PaperPositionMode','auto')    
    print('-painters','-r300','-dpng',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_orient' int2str(sliceThisDim) '_slice' int2str(kk)])
end

%% Subject 1 Signals from voxels in ROIs

% plot the slice with ROIs
curPos = [-1,1,-20]; 
sliceThisDim = 1; 
imDims = [-90 -120 -120; 90 130 90];
overlayPlot = ppgR;
cod_th = 0.6;
bbOverlayFuncAnat(overlayPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,1,1,cod_th);

bb_roi = bb_subs_rois(ss); % script that lists voxel indices for subject 1

% Loop through ROI indices/labels
for roi_ind = 1:8
    curPos = bb_roi(roi_ind).curPos;
    voxelLabel = bb_roi(roi_ind).voxelLabel;
    plot(curPos(2),curPos(3),'w.')
%     text(curPos(2),curPos(3),voxelLabel,'Color',[1 1 1])
end
% print('-painters','-r300','-dpng',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_orient' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim)) '_ROIs'])
% print('-painters','-r300','-dpng',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_orient' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim)) '_ROItext'])
%%

% load COD
ppgRname = [save_name_base '_codPPG.nii.gz'];

% ppg locked time series
ppgTS = niftiRead([save_name_base '_PPGtrigResponse.nii.gz']); 
% associated time t
load([save_name_base '_PPGtrigResponseT.mat'],'t');
% load std of time series
ppgTSstd = niftiRead([save_name_base '_PPGtrigResponse_std.nii.gz']); % ppg triggered time series

bb_roi = bb_subs_rois(ss); % script that lists voxel indices

for roi_ind = 1:8
    % Call ROI indices/labels
    curPos = bb_roi(roi_ind).curPos;
    voxelLabel = bb_roi(roi_ind).voxelLabel;
    % get the timeseries
    [voxelTs] = bbGetVoxelTimeseries(ppgTS,acpcXform,curPos);
    [voxelTsStd] = bbGetVoxelTimeseries(ppgTSstd,acpcXform,curPos);

    % plot the timeseries
%     figure('Position',[0 0 150 150]),hold on
    figure('Position',[0 0 150 200]),hold on
    plot(t,zeros(size(t)),'Color',[.5 .5 .5])
    plot([0 0],[-2 2],'Color',[.5 .5 .5])
    upErr = 100*squeeze(voxelTs) + 100*voxelTsStd; % mean + standard deviation
    lowErr = 100*squeeze(voxelTs) - 100*voxelTsStd;
    fill([t t(end:-1:1)],[upErr; lowErr(end:-1:1)],[.5 .5 .5],'EdgeColor',[.5 .5 .5])
    plot(t,100*squeeze(voxelTs),'k','LineWidth',2)
    axis tight
    ylabel('% signal modulation') % (signal - mean)./mean
    xlabel('time (s)')
    title(bb_roi(roi_ind).voxelLabel)
    set(gcf,'PaperPositionMode','auto')
%     print('-painters','-r300','-dpng',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_TraceOnAnat_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3)) '_' voxelLabel])
%     print('-painters','-r300','-depsc',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_TraceOnAnat_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3)) '_' voxelLabel])
end

%%
%%
%% %% Overlay 2: timeseries and anatomy
%%
%%


ppgTSplot = ppgTS;
ppgTSplot.data(:,:,:,t<-.1 | t>1)=[]; % plot these times from curve

% Scale time series amplitude by R, plots are generated with respect to the maximum.
niColor = ppgR; % use R2 map to scale later: r2_scale=sqrt(imgSlice2.^2);
maxTS = max(abs(ppgTSplot.data),[],4); % get the max of each curve
ppgTSplot.data = bsxfun(@rdivide,ppgTSplot.data,maxTS); % devide by the max of each curve (sets all curves to 1 max)
ppgTSplot.data = bsxfun(@times,ppgTSplot.data,niColor.data); % multiply by r^2 to set less reliable curves to zero

bbOverlayTimeseriesAnat(ppgTSplot,niColor,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)

clear niColor ppgTSplot


%% plot MRV and timeseries
% Get the MRVenogram:
niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));
% load coregistration matrix (for the venogram):
xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));

bbOverlayTimeseriesVeno(ppgTSplot,niColor,niVeno,acpcXform,xf_veno.acpcXform,sliceThisDim,imDims,curPos)

