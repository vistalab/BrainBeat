clear all
close all

% script produced transformation matrix between the anatomy and functional
% scans, walk through step-by-step

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% chdir(dDir)

%%
%% Overlay anat and func
%%

%% The T2* data are here.  

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

in_data = 'PPG';
s_nr = 1;
scan_nr = 1;

s_info = bb_subs(s_nr);
subj = s_info.subj;
    
% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

% % Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));

scan = s_info.scan{scan_nr};
scanName = s_info.scanName{scan_nr};

% loaf the fMRI data
% ni = niftiRead(fullfile(dDir,subj,scan, [scanName '.nii.gz']));

% load coregistration matrix (for the functionals):
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

%%%% Overlay 1: cod and anatomy
ppgRname = fullfile(dDir,subj,scan,[scanName '_cod' in_data '.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

%% Subject 1 Sagittal slices T1 + COD
curPos = [-10,1,-20]; 
sliceThisDim = 1; 
imDims=[-90 -120 -120; 90 130 90];
overlayPlot = ppgR;
cod_th = 0.6;
for kk = -60:10:60
    curPos(1) = kk;
    bbOverlayFuncAnat(overlayPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,1,1,cod_th);
    title(['slice ' int2str(kk) ', R^2>' num2str(cod_th,3)])
    set(gcf,'PaperPositionMode','auto')    
    print('-painters','-r300','-dpng',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_orient' int2str(sliceThisDim) '_slice' int2str(kk)])
end

%% Subject 1 Axial slices T1 + COD
curPos = [-10,1,-20]; 
sliceThisDim = 3; 
imDims=[-90 -120 -120; 90 130 90];
overlayPlot = ppgR;
cod_th = 0.6;
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

bb_roi = bb_subs_rois(s_nr); % script that lists voxel indices for subject 1

% Loop through ROI indices/labels
for roi_ind = 1:8
    curPos = bb_roi(roi_ind).curPos;
    voxelLabel = bb_roi(roi_ind).voxelLabel;
    plot(curPos(2),curPos(3),'w.')
    text(curPos(2),curPos(3),voxelLabel,'Color',[1 1 1])
end
% print('-painters','-r300','-dpng',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_orient' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim)) '_ROIs'])
print('-painters','-r300','-dpng',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_orient' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim)) '_ROItext'])
%%
% load time series and associated time
ppgTSname = fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponseT.mat']),'t');

% load error of time series and associated time
ppgTSnameError = fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponse_std.nii.gz']);
ppgTSError = niftiRead(ppgTSnameError); % ppg triggered time series

bb_roi = bb_subs_rois(s_nr); % script that lists voxel indices

for roi_ind = 1:8
    % Call ROI indices/labels
    curPos = bb_roi(roi_ind).curPos;
    voxelLabel = bb_roi(roi_ind).voxelLabel;
    % get the timeseries
    [voxelTs] = bbGetVoxelTimeseries(ppgTS,acpcXform,curPos);
    [voxelTsStd] = bbGetVoxelTimeseries(ppgTSError,acpcXform,curPos);

    % plot the timeseries
    figure('Position',[0 0 150 150]),hold on
    plot(t,zeros(size(t)),'Color',[.5 .5 .5])
    plot([0 0],[-2 2],'Color',[.5 .5 .5])
    upErr = 100*squeeze(voxelTs) + 100*voxelTsStd; % mean + 2 standard error
    lowErr = 100*squeeze(voxelTs) - 100*voxelTsStd;
    fill([t t(end:-1:1)],[upErr; lowErr(end:-1:1)],[.5 .5 .5],'EdgeColor',[.5 .5 .5])
    plot(t,100*squeeze(voxelTs),'k','LineWidth',2)
    axis tight
    ylabel('% signal modulation') % (signal - mean)./mean
    xlabel('time (s)')
    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_TraceOnAnat_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3)) '_' voxelLabel])
    print('-painters','-r300','-depsc',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_TraceOnAnat_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3)) '_' voxelLabel])
end

%%
%%
%% %% Overlay 2: timeseries and anatomy
%%
%%

% load time series and associated time
ppgTSname = fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponseT.mat']),'t');

% load correlation between even and odd scans for colors of timeseries

ppgTSplot = ppgTS;
if isequal(in_data,'PPG')
    ppgTSplot.data(:,:,:,t<-.1 | t>1)=[]; % plot these times from curve
elseif isequal(in_data,'RESP')
    ppgTSplot.data(:,:,:,t<-.2 | t>4)=[]; % plot these times from curve
end
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

