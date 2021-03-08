

%% The functionals were not writted in the same space as the T1 
% here save the coregistration matrix to the f.._codPPG.nii in sto_xyz/ijk
% this is necessary for the normalization

% dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

in_data = 'PPG';
s = 1;
scan_nr = 3;

s_info = bb_subs(s);
subj = s_info.subj;
    
scan = s_info.scan{scan_nr};
scanName = s_info.scanName{scan_nr};

fcod = fullfile(dDir,subj,scan, ['f' scanName '_codPPG.nii']);
ni = niftiRead(fcod);

% load coregistration matrix for the functionals
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% save this coregistration matrix in the f_...codPPG.nii
ni.qto_xyz = acpcXform;
ni.qto_ijk = inv(acpcXform);
ni.sto_xyz = acpcXform;
ni.sto_ijk = inv(acpcXform);

niftiWrite(ni,fcod);


%% now we can normalize the cod image

spm('Defaults','fmri')

% A) volume where parameters are estimated: anat
% code will search for y_anat.nii file to use for normalization
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,['f' s_info.anatName '.nii']));

% B) volume you want to normalize, needs to be coregistered with A)
% code will append a w before this file
fcod = fullfile(dDir,subj,scan, ['f' scanName '_codPPG.nii']);

flags.preserve  = 0;
flags.bb        = [-90 -120 -60; 90 96 130];
flags.vox       = [1 1 1]; % here is the voxel size
flags.interp    = 0;
flags.wrap      = [0 0 0];
flags.prefix    = 'w';

job.subj.vol{1} = niAnatomy;
job.subj.resample{1} = fcod;
job.woptions = flags;

% normalize the image with electrodess
spm_run_norm(job);


