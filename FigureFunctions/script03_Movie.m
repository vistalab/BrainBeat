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

%% Movie with overlay
movieName = [dDir 'movies/PPGtimeseries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))];

bbOverlayFuncAnatMovie(ppgTSplot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,movieName,t)

%% Movie with dots
movieName = [dDir 'movies/PPGtimeseries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim)) '_dots'];

bbOverlayDotsAnatMovie(ppgTSplot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,movieName,t)

%% Movie with SVD

data_in = 'PPG';

% load average of even and odd heartbeats:
ppgTSodd=niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_odd.nii.gz']));
ppgTSeven=niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_even.nii.gz']));

% Set maximum of ppgTS to 1 for each voxel
ppgTSodd.data = ppgTSodd.data ./ repmat(max(abs(ppgTSodd.data),[],4),[1,1,1,size(ppgTSodd.data,4)]);
ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
% Multiply by correlation size (absolute)
ppgTSodd.data = ppgTSodd.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSodd.data,4)]));
ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));

% reshape to voxel X time
a = reshape(ppgTSodd.data,[numel(ppgTSodd.data(:,:,:,1)) length(t)]);

% select times to include in SVD
if isequal(data_in,'PPG')
    a = a(:,t>=-.5 & t<=2);
    t_sel = t(t>=-.5 & t<=2);
elseif isequal(data_in,'RESP')
    a = a(:,t>=-.5 & t<=4);
    t_sel = t(t>=-.5 & t<=4);
end
meanTS = mean(a,2);
a = a-repmat(meanTS,1,size(a,2)); % subtract the mean
[u,s,v]=svd(a','econ');
s=diag(s);

%%%% test for the sign of the 2nd component (also check for 1st???)
% in the 2nd component, the first peak should be negative
[~,pm_i]=findpeaks(double(u(:,2)),'minpeakdistance',10);
[~,nm_i]=findpeaks(-double(u(:,2)),'minpeakdistance',10);
first_peak = min([pm_i; nm_i]);
if u(first_peak,2)<0
    % do nothing
elseif u(first_peak,2)>0
    % reverse sign pc and weights
    u(:,2)=-u(:,2); v(:,2) = -v(:,2);
end

% get to cumulative explained variance:
var_explained=cumsum(s.^2) / sum(s.^2);

out(1).weights = NaN(size(ppgTSodd.data,1),size(ppgTSodd.data,2),size(ppgTSodd.data,3));
out(2).weights = NaN(size(ppgTSodd.data,1),size(ppgTSodd.data,2),size(ppgTSodd.data,3));
for k=1:2
    out(k).weights = reshape(v(:,k),[size(ppgTSodd.data,1) size(ppgTSodd.data,2) size(ppgTSodd.data,3)]);
end

%%%%% MODEL WITH 2 COMPONENTS:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';
svdResults.model = reshape(pred,[size(ppgTSodd.data,1) size(ppgTSodd.data,2) size(ppgTSodd.data,3) size(ppgTSodd.data,4)]);

%%%%% MODEL ERROR
% train and test-sets
train_set = reshape(ppgTSodd.data,[prod(ppgTSodd.dim(1:3)) ppgTSodd.dim(4)]);
test_set = reshape(ppgTSeven.data,[prod(ppgTSeven.dim(1:3)) ppgTSeven.dim(4)]);
% test-retest error
test_train_error = sqrt(sum((test_set - train_set).^2,2));
% model error
test_model_error = sqrt(sum((test_set - pred).^2,2));
% relative RMS error:
rel_rms_error = test_model_error./test_train_error;

svdResults.error = reshape(rel_rms_error,[size(ppgTSodd.data,1) size(ppgTSodd.data,2) size(ppgTSodd.data,3)]);

%%
physio  = physioCreate('nifti',ni);
ppgRate = physioGet(physio,'ppg rate'); % peaks per second

tMin = find(t>=-.5,1); % first sample after 1 heartbeat
tMax = find(t>1./ppgRate-.5,1); % first sample after 1 heartbeat

tSelect = tMin:tMax;

clear tMin tMax
%% Negative values - veins and arteries:

niPlot = ppgTSplot;
% get the times for 1 heartbeat cycle, we can then loop the movie later
tPlot = t(tSelect);
niPlot.data = niPlot.data(:,:,:,tSelect);

% Set COD threshold
for kk = 1:size(niPlot.data,4)
    currentVol = niPlot.data(:,:,:,kk);
    currentVol(ppgR.data<.6) = 0;
    niPlot.data(:,:,:,kk) = currentVol;
    clear currentVol
end

% Only make movie of positive or negative PC1
for kk = 1:size(niPlot.data,4)
    currentVol = niPlot.data(:,:,:,kk);
    currentVol(out(1).weights>0) = 0;
    niPlot.data(:,:,:,kk) = currentVol;
    clear currentVol
end
movieName = [dDir 'movies/PPGtimeseries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim)) '_neg2'];
bbOverlayDotsAnatMovie(niPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,movieName,tPlot)

% Positive values - CSF, plot with peak in same direction as negative ones, b/c of colorscale:

niPlot = ppgTSplot;
% get the times for 1 heartbeat cycle, we can then loop the movie later
tPlot = t(tSelect);
niPlot.data = -niPlot.data(:,:,:,tSelect);

% Set COD threshold
for kk = 1:size(niPlot.data,4)
    currentVol = niPlot.data(:,:,:,kk);
    currentVol(ppgR.data<.6) = 0;
    niPlot.data(:,:,:,kk) = currentVol;
    clear currentVol
end

% Only make movie of positive or negative PC1
for kk = 1:size(niPlot.data,4)
    currentVol = niPlot.data(:,:,:,kk);
    currentVol(out(1).weights<0) = 0;
    niPlot.data(:,:,:,kk) = currentVol;
    clear currentVol
end
movieName = [dDir 'movies/PPGtimeseries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim)) '_pos2'];
bbOverlayDotsAnatMovie(niPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,movieName,tPlot)

