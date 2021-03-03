clear all
close all

% script produced transformation matrix between the anatomy and functional
% scans, walk through step-by-step

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
% dDir = '/biac4/wandell/data/BrainBeat/data';
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

%%
curPos = [-5,1,-30]; 
sliceThisDim = 3; 
% imDims=[-90 -120 -60; 90 130 90];
imDims=[-90 -120 -120; 90 130 90];

overlayPlot = ppgR;

bbOverlayFuncAnat(overlayPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,1,1,.2);


%%
%%%% Overlay 2: timeseries and anatomy
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


%%
%%
%% The rest seems to be unnecessary?
%%
%%
%%


%%
%% Allign functional to a good functional with SPM
%%
% clear all
% close all
% dDir = '/biac4/wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

s_nr = 4;
scan_nr = 9;
ref_scan_nr = 3;

s_info = bb_subs(s_nr);
subj=s_info.subj;

% Get the new ref scan:
%     niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.scan{ref_scan_nr}, [s_info.scanName{ref_scan_nr} '.nii.gz']));
niAnatomy.data = niAnatomy.data(:,:,:,1);

%%%%% coregister the functionals to the ref funx:

scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};

fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
ni = niftiRead(fmri);

%%%%% use the first nifti to align, this one has the most structural info:
ni1=ni;
ni1.data=ni1.data(:,:,:,1);
ni1.dim=ni1.dim(1:3);
ni1.pixdim=ni1.pixdim(1:3);
niAnatomy.pixdim=niAnatomy.pixdim(1:3); % only uses the first three dimensions

% align functionals to the T1:
acpcXform = dtiRawAlignToT1(ni1,niAnatomy,[], [], false, 1); % last 1 adds a figure
% this saves the realignment matrix in the folder of the functionals, 

% load ref scan acpc x-form
ref_acpc = load(fullfile(dDir,subj,s_info.scan{ref_scan_nr},...
    [s_info.scanName{ref_scan_nr} 'AcpcXform_new.mat']));

% now fix the acpcXform such that it goes to T1 space
acpcXform_new = ref_acpc.acpcXform_new * niAnatomy.qto_ijk * acpcXform; % funx -> ref xyz -> ref ijk -> ref acpc
save(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']),'acpcXform_new')



