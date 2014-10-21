%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%
%

%% Base data directory on a Mac mounting biac4 (wandell's machine)
dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';

% chdir(dDir)

%% The T2* data are here.  

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

subj ='20141017_1242';    % Date _ Time out of NIMS
scan ='6_1_mux8fov4_r1_25s_4mmFA25';  % A particular data set
fmri = fullfile(dDir,subj,scan,'8202_6_1.nii.gz');
ni = niftiRead(fmri);


%% We should get the anatomicals up too at some point

anat      ='9_1_T1w_1mm_sag';   % Anatomical data
anat      = fullfile(dDir,subj,anat,'8202_9_1.nii.gz');
niAnatomy = niftiRead(anat);

showMontage(niAnatomy.data);

% mrViewer should be able to take a nifti structure as input data
mrViewer(anat)

%% Let's pick a physiology file and do something


%% End