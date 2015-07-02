clear all
close all

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
dDir = '/biac4/wandell/data/BrainBeat/data';

% chdir(dDir)

%% The T2* data are here.  

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

s_nr = 3;
s_info = bb_subs(s_nr);
subj=s_info.subj;


%% Get the anatomicals:

niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

%% Get the MRVenogram:

niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii.gz']));

%%
%% coregister the venogram to the T1:
%%

% if already done, load rotation matrix:
xf_veno=load(fullfile(dDir,subj,veno,[veno_code 'AcpcXform.mat']));

%%
%% now do an SVD on the responses:
%%

% load functional data
scan_nr = 3;
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};
fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
ni = niftiRead(fmri);
% load coregistration matrix:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform.mat']))

% Load the correlation with heartbeat (made with bbCorrelate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_corrPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

% Load the timeseries around PPG peak (made with bbResponse2physio):
ppgTSname = fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT.mat']),'t');

a = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);
a = a(:,t>=-0.5 & t<=2);
t_sel = t(t>=-0.5 & t<=2);
% a = a-repmat(mean(a,2),1,size(a,2)); % subtract the mean
[u,s,v]=svd(a','econ');
s=diag(s);

% get to cumulative explained variance:
var_explained=cumsum(s.^2) / sum(s.^2)

%%

phy = bbGet(ni,'physio');
ppgRate = physioGet(phy,'PPGrate');

t_in = t(t>=-0.5 & t<=.5);
ts_in = double(a(436,t>=-0.5 & t<=.5));

my_options = optimset('Display','off','Algorithm','trust-region-reflective');
[x] = lsqnonlin(@(x) fit_func_Gaussians(x,ts_in,t_in),...
    [0 0 .1 .5 .25 .1 -1],[0 t_in(1) 0 0 t_in(1) 0 -Inf],[Inf t_in(end) 1 Inf t_in(end) 1 Inf],... %X0,LB,UB
    my_options);


fitted_line = x(7) + ...
    x(1)*x(3)*sqrt(2*pi)*normpdf(t_in,x(2),x(3)) + ...
    x(4)*x(6)*sqrt(2*pi)*normpdf(t_in,x(5),x(6));

figure,hold on
plot(t_in,ts_in,'r')
plot(t_in,fitted_line,'k')

%%

