s_nr = 1;

% Load PPG responses:
s_info      = bb_subs(s_nr);
subj        = s_info.subj;
scan        = s_info.scan{scan_nr};
scanName    = s_info.scanName{scan_nr};

% % Get the anatomicals:
% niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

t1_name = fullfile(dDir,['sub-' int2str(s_nr)],'anat',[s_info.anatName '.nii.gz']);
t1_save = fullfile(dDir,['sub-' int2str(s_nr)],'anat',['sub-' int2str(s_nr) '_T1w.nii.gz']);

mrAnatAverageAcpcNifti(t1_name, t1_save);
