clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';


%% The T2* data are here.  

% Select a subject and scan nummer
s_nr = 2;
scan_nr = 3;

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

anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii']);
niAnatomy = niftiRead(anat);

%% Overlay between functionals and anatomy

sliceThisDim = 1;
if s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-1,10,-20];
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [1,4,38];
end

movieName = [dDir 'movies/PPGtimeseries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))];

load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT']),'t')

% load average of all odd heartbeats:
ppgTS = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse.nii.gz']));
ppgTSplot = ppgTS;
ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG
niColor = ppgR; % use R2 map to scale later: r2_scale=sqrt(imgSlice2.^2);

% Scale time series amplitude by R, plots are generated with respect to the maximum.
maxTS = max(abs(ppgTSplot.data),[],4); % get the max of each curve
ppgTSplot.data = bsxfun(@rdivide,ppgTSplot.data,maxTS); % devide by the max of each curve (sets all curves to 1 max)
ppgTSplot.data = bsxfun(@times,ppgTSplot.data,niColor.data); % multiply by R to set less reliable curves to zero

% Note that this only shows positive values
bbOverlayFuncAnatMovie(ppgTSplot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,movieName,t)

%% 
movieName = [dDir 'movies/PPGtimeseries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim)) '_dots'];

bbOverlayDotsAnatMovie(ppgTSplot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,movieName,t)

% add a threshold for dots to plot
