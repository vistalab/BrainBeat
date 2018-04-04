clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';
% chdir(dDir)

%% The T2* data are here.  

% Select a subject and scan nummer
s_nr = 4;
scan_nr = 9;

subs = bb_subs(s_nr);
subj = subs.subj;
scan = subs.scan{scan_nr};
scanName = subs.scanName{scan_nr};
fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
if ~exist(fmri,'file')
    clear ni
    error('filename %s does not exist',fmri)
end
ni = niftiRead(fmri);

% load coregistration matrix for the functionals:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

% Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,subs.veno,[subs.venoName '.nii']));
% % load coregistration matrix (for the venogram):
% xf_veno=load(fullfile(dDir,subj,subs.veno,[subs.venoName 'AcpcXform.mat']));

%% Anatomicals

anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii.gz']);
niAnatomy = niftiRead(anat);

%% quick overlay between functionals and anatomy

sliceThisDim = 3;
if s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [1,10,-20];
%     curPos = [-29,-14,-50]; % 03
%     curPos = [-3,30,-43]; % 03
    curPos = [1,10,-20];
    curPos = [1,1,-20];
    curPos = [-11 34 -71]; % Carotid
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
%     curPos = [0,4,38];
elseif s_nr == 4
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
end
niFunc = ni;

% niFunc.data = ni.data(:,:,:,1); % overlay the first functional - more structure visible
% bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
% title('First functional on anatomy')
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/checkCoreg/' subj '_' scan '_Func1onAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

niFunc.data = mean(ni.data(:,:,:,5:end),4); % overlay the mean functional
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
title('Mean functional on anatomy')
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/checkCoreg/' subj '_' scan '_MeanFuncOnAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])


%% Plot physiology data

figure('Position',[0 0 600 600])

% get data
physio      = physioCreate('nifti',ni);
ppgData     = physioGet(physio,'ppg data');
respData    = physioGet(physio,'resp data');

% plot PPG peaks on PPG signal
subplot(3,4,1:3)
ppgPeaks    = physioGet(physio,'ppg peaks');
tPPG = physioGet(physio,'ppg sample times');
srate = physioGet(physio,'ppg srate');
ppgPeakSamples = round(ppgPeaks*srate);
plot(ppgPeaks,ppgData(ppgPeakSamples),'ro',tPPG,ppgData,'k-');
xlabel('time (s)'); grid on
title('PPG')
xlim([0 10])

% plot one PPG curve
subplot(3,4,4)
ppgCurve = physioGet(physio,'ppg ppgcurve');
ppgCurveT = physioGet(physio,'ppg ppgtcurve');
ppgSrate = physioGet(physio,'PPGsrate'); %  = physio.ppg.srate
plot(ppgCurveT,ppgCurve,'k','LineWidth',2);
xlim([ppgCurveT(1) ppgCurveT(end)])
xlabel('time (s)')
% get PPG rate
ppgRate = physioGet(physio,'ppg rate');
title(['PPG rate = ' num2str(ppgRate)])

% plot RESP peaks on RESP signal
subplot(3,4,5:7)
respPeaks = physioGet(physio,'resp peaks');
tResp = physioGet(physio,'resp sample times');
srate = physioGet(physio,'resp srate');
respPeakSamples = round(respPeaks*srate);
plot(respPeaks,respData(respPeakSamples),'ro',tResp,respData,'k-');
xlabel('time (s)'); grid on
xlim([0 50])
title('RESP')

% plot one RESP curve
subplot(3,4,8)
respCurve = physioGet(physio,'resp respcurve');
respCurveT = physioGet(physio,'resp resptcurve');
respSrate = physioGet(physio,'RESPsrate');
plot(respCurveT,respCurve,'k','LineWidth',2);
xlabel('time (s)')
xlim([respCurveT(1) respCurveT(end)])
% get RESP rate
respRate = physioGet(physio,'resp rate');
title(['RESP rate = ' num2str(respRate)])

% plot PPG signal with scan onsets
subplot(3,1,3),hold on
plot(ppgPeaks,ppgData(ppgPeakSamples),'ro',tPPG,ppgData,'k-');
plot(tPPG(physio.ppg.scan_onset==1),0,'b*')
title('PPG signal (black) and scan onsets (blue)')
xlim([0 10])

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-depsc',[dDir './figures/physio/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_physioTrace'])
% print('-painters','-r300','-dpng',[dDir './figures/physio/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_physioTrace'])


%% Plot MRI data

in_data = 'PPG';

sliceThisDim = 1;
if s_nr == 1 % subject number
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-1 26 -63]; 
elseif s_nr == 2 % subject number
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [-10 50 -21]; % left lateral ventricle - frontal site
    curPos = [-1 72 -19]; % ExampleSite ACA
elseif s_nr == 3 %subject number
    imDims = [-90 -120 -100; 90 130 110];
%     curPos = [1,4,38];
%     curPos = [-14 24 -16]; % LCarotid
    curPos = [1 11 -16]; % Basilar
elseif s_nr == 4 %subject number
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [-5,4,38];
end

% load time series and associated time
ppgTSname = fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponseT.mat']),'t');

% load correlation between even and odd scans for colors of timeseries
% ppgRname = fullfile(dDir,subj,scan,[scanName '_corr' in_data '.nii.gz']);
% ppgR = niftiRead(ppgRname); % correlation with PPG
% ppgR.data = ppgR.data.^2;
ppgRname = fullfile(dDir,subj,scan,[scanName '_cod' in_data '.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

% %%%% Overlay 1: functionals and anatomy
% niFunc = ni;
% niFunc.data = mean(ni.data(:,:,:,5:end),4); % overlay the mean functional
% bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)

%%%% Overlay 2: timeseries and anatomy
ppgTSplot = ppgTS;
if isequal(in_data,'PPG')
    ppgTSplot.data(:,:,:,t<-.1 | t>1)=[]; % plot these times from curve
elseif isequal(in_data,'RESP')
    ppgTSplot.data(:,:,:,t<-.2 | t>4)=[]; % plot these times from curve
end

% Scale factor for size and color:
niColor = ppgR; % use R2 map to scale later: r2_scale=sqrt(imgSlice2.^2);

% Scale time series amplitude by R, plots are generated with respect to the maximum.
maxTS = max(abs(ppgTSplot.data),[],4); % get the max of each curve
ppgTSplot.data = bsxfun(@rdivide,ppgTSplot.data,maxTS); % devide by the max of each curve (sets all curves to 1 max)
ppgTSplot.data = bsxfun(@times,ppgTSplot.data,niColor.data); % multiply by r^2 to set less reliable curves to zero

bbOverlayTimeseriesAnat(ppgTSplot,niColor,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
% the plotted X and Y do not correspond to actual xyz/mm coordinates, but
% are in a different frame
title('Colors and size scaled by the COD(R)')

clear niColor ppgTSplot
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/voxelTimeSeries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_TraceOnAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

% bbOverlayTimeseriesVeno(ppgTSplot,niColor,niVeno,acpcXform,xf_veno.acpcXform,sliceThisDim,imDims,curPos)

%%
figure('Position',[0 0 100,200]),hold on
cm = colormap(jet);
for kk = 1:size(cm,1)
    plot(1,kk,'.','Color',cm(kk,:),'MarkerSize',20)
end
axis off
set(gcf,'PaperPositionMode','auto')
% 
% print('-painters','-r300','-dpng',[dDir './figures/voxelTimeSeries/TScolorbar'])
% print('-painters','-r300','-depsc',[dDir './figures/voxelTimeSeries/TScolorbar'])

%% Get voxel timeseries

% Now say we have a voxel let's pull out the raw timeseries and show
% the processing we did:

in_data = 'PPG';

bb_roi = bb_subs_rois(s_nr); % script that lists voxel indices

for roi_ind = 1%1:length(bb_roi)

if s_nr == 2 % subject number
    imDims = [-90 -120 -120; 90 130 90];
       
    % Enter Cursor Position for Voxel
%     curPos = [-1 21 -86]; 
%     voxelLabel = '';

    curPos = bb_roi(roi_ind).curPos;
    voxelLabel = bb_roi(roi_ind).voxelLabel;

%     curPos = [-11 34 -71]; 
%     voxelLabel = 'LCarotid1';
%     curPos = [-11 36 -61]; 
%     voxelLabel = 'LCarotid2';
%     curPos = [8 36 -60]; 
%     voxelLabel = 'RCarotid1';
elseif s_nr == 3 %subject number
    imDims = [-90 -120 -100; 90 130 110];
%     curPos = [1,4,38];
    curPos = [-14 24 -16]; % LCarotid
    voxelLabel = 'LCarotid1';
%     curPos = [-1 11 -16]; % Basilar
%     voxelLabel = 'Basilar1';
%     curPos = bb_roi(roi_ind).curPos;
%     voxelLabel = bb_roi(roi_ind).voxelLabel;
elseif s_nr == 4 %subject number
    imDims = [-90 -120 -100; 90 130 110];
%     curPos = [1,4,38];
%     curPos = [-14 24 -16]; % LCarotid
%     voxelLabel = 'LCarotid1';
    curPos = bb_roi(roi_ind).curPos;
    voxelLabel = bb_roi(roi_ind).voxelLabel;

end

% load time series and associated time
ppgTSname = fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponseT.mat']),'t');

% load error of time series and associated time
ppgTSnameError = fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponse_std.nii.gz']);
ppgTSError = niftiRead(ppgTSnameError); % ppg triggered time series

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
print('-painters','-r300','-dpng',[dDir './figures/voxelTimeSeries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_TraceOnAnat_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3)) '_' voxelLabel])
print('-painters','-r300','-depsc',[dDir './figures/voxelTimeSeries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_TraceOnAnat_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3)) '_' voxelLabel])

end

%% plot Raw timeseries with cardiac signal

curPos = [-10 50 -21]; % subj = 2, scan_nr = 3
% curPos = [-14 24 -16];

% get the timeseries 
[voxelTs] = bbGetVoxelTimeseries(ni,acpcXform,curPos);

% plot the timeseries
figure('Position',[0 0 300 150]),hold on
plot(t,zeros(size(t)),'Color',[.5 .5 .5])
plot([0 0],[-2 2],'Color',[.5 .5 .5])
total_t = [1:length(voxelTs)]./(1./ni.pixdim(4));
plot(total_t,voxelTs,'k','LineWidth',1)
ylabel('raw fMRI TS') % (signal - mean)./mean
xlabel('time (s)')
xlim([total_t(1) total_t(end)])
set(gcf,'PaperPositionMode','auto')

physio      = physioCreate('nifti',ni);
tPPG        = physioGet(physio,'ppg sample times');
ppgNorm     = zscore(physio.ppg.data);
ppgNorm     = ppgNorm + mean(voxelTs);
plot(tPPG,ppgNorm,'r')

% xlim([24 29])
% ylim([63 73])
% print('-painters','-r300','-dpng',[dDir './figures/voxelTimeSeries/' subj '_' scan '_RawTrace_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3))])
% print('-painters','-r300','-depsc',[dDir './figures/voxelTimeSeries/' subj '_' scan '_RawTrace_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3))])

