clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%

%% Base data directory on a Mac mounting biac4 (wandell's machine)

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';
% chdir(dDir)

%% The T2* data are here.  

% Select a subject and scan nummer
s_nr = 2;
scan_nr = 2;

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

%% Anatomicals

anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii']);
niAnatomy = niftiRead(anat);

%% quick overlay between functionals and anatomy
sliceThisDim = 1;
if s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [1,10,-20];
%     curPos = [-29,-14,-50]; % 03
%     curPos = [-3,30,-43]; % 03
    curPos = [1,10,-20];
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
%     curPos = [0,4,38];
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
print('-painters','-r300','-dpng',[dDir './figures/checkCoreg/' subj '_' scan '_MeanFuncOnAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])


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
ppgCurve = physioGet(physio,'ppg curve');
ppgSrate = physioGet(physio,['PPGsrate']); %  = physio.ppg.srate
plot((1:length(ppgCurve))/ppgSrate,ppgCurve,'k','LineWidth',2);
xlim([0 max((1:length(ppgCurve))/ppgSrate)])
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
respCurve = physioGet(physio,'resp curve');
respSrate = physioGet(physio,['RESPsrate']);
plot((1:length(respCurve))/respSrate,respCurve,'k','LineWidth',2);
xlabel('time (s)')
xlim([0 max((1:length(respCurve))/respSrate)])
% get RESP rate
respRate = physioGet(physio,'resp rate');
title(['RESP rate = ' num2str(respRate)])

% plot PPG signal with scan onsets
subplot(3,1,3),hold on
plot(ppgPeaks,ppgData(ppgPeakSamples),'ro',tPPG,ppgData,'k-');
plot(tPPG(physio.ppg.scan_onset==1),0,'b*')
title('PPG signal (black) and scan onsets (blue)')
xlim([0 10])

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-depsc',[dDir './figures/physio/' subj '_' scan '_physioTrace'])
% print('-painters','-r300','-dpng',[dDir './figures/physio/' subj '_' scan '_physioTrace'])

%% Plot MRI data

in_data = 'PPG';

sliceThisDim = 1;
if s_nr == 2 % subject number
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-10 50 -21]; % left lateral ventricle - frontal site
%     curPos = [-3,40,-48]; % optic chiasm
%     curPos = [21,10,-80]; % right trigeminal?
elseif s_nr == 3 %subject number
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [1,4,38];
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
% print('-painters','-r300','-dpng',[dDir './figures/checkCoreg/' subj '_' scan '_TraceOnAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])
% print('-painters','-r300','-dpng',[dDir './figures/voxelTimeSeries/' subj '_' scan '_TraceOnAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

%%
figure('Position',[0 0 100,200]),hold on
cm = colormap(jet);
for kk = 1:size(cm,1)
    plot(1,kk,'.','Color',cm(kk,:),'MarkerSize',20)
end
axis off
set(gcf,'PaperPositionMode','auto')

print('-painters','-r300','-dpng',[dDir './figures/voxelTimeSeries/TScolorbar'])
print('-painters','-r300','-depsc',[dDir './figures/voxelTimeSeries/TScolorbar'])

%% Get voxel timeseries

% Now say we have a voxel let's pull out the raw timeseries and show
% the processing we did:

in_data = 'PPG';

% sliceThisDim = 1;
if s_nr == 2 % subject number
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [-10 50 -21]; % left lat ventricle - in mm coordinates on the anatomical MRI
    curPos = [-1,40,-48]; % 
elseif s_nr == 3 %subject number
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [1,4,38];
end

% load time series and associated time
ppgTSname = fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponseT.mat']),'t');

% get the timeseries
[voxelTs] = bbGetVoxelTimeseries(ppgTS,acpcXform,curPos);

% plot the timeseries
figure('Position',[0 0 150 150]),hold on
plot(t,zeros(size(t)),'Color',[.5 .5 .5])
plot([0 0],[-2 2],'Color',[.5 .5 .5])
plot(t,100*squeeze(voxelTs),'k','LineWidth',2)
xlim([t(1) t(end)])
ylabel('% signal modulation') % (signal - mean)./mean
xlabel('time (s)')
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/voxelTimeSeries/' subj '_' scan '_TraceOnAnat_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3))])
% print('-painters','-r300','-depsc',[dDir './figures/voxelTimeSeries/' subj '_' scan '_TraceOnAnat_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3))])

%% plot Raw timeseries with cardiac signal

curPos = [-10 50 -21]; % subj = 2, scan_nr = 3

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

xlim([24 29])
ylim([63 73])
print('-painters','-r300','-dpng',[dDir './figures/voxelTimeSeries/' subj '_' scan '_RawTrace_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3))])
print('-painters','-r300','-depsc',[dDir './figures/voxelTimeSeries/' subj '_' scan '_RawTrace_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3))])

