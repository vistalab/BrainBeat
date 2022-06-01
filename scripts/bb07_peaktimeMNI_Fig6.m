
% This script generates Figure 6 from the manuscript titled:
%
% Measuring brain beats: cardiac-aligned fast fMRI signals
% Dora Hermes, Hua Wu, Adam B. Kerr, Brian Wandell
%


clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%% Normalize PC1/PC2 and the slope, latency and width to MNI152

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,1,1,1,1};

for ss = 1:length(sub_labels)
    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);
    
    %%%%%% now we can normalize the PC1 and PC2 beta weight image
    spm('Defaults','fmri')
    
    % A) volume where parameters are estimated: anat
    % code will search for y_anat.nii file to use for normalization
    niAnatomy = fullfile(dDir,'derivatives','spmSegmentation',['sub-' sub_label],['ses-' ses_label],...
                ['sub-' sub_label '_ses-' ses_label '_T1w']);
    
    % B) volume you want to normalize, needs to be coregistered with A)
    % code will append a w before this file
    pc1_NameSave = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_canoPc1Weights.nii']);
    pc2_NameSave = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_canoPc2Weights.nii']);
    
    % Get base name
    slope_NameSave = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_modelSlope.nii']);
    onset_NameSave = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_modelOnset.nii']);
    fwhm_NameSave = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_model_FWHM.nii']);
    
    flags.preserve  = 0;
    flags.bb        = [-90 -120 -60; 90 96 130];
    flags.vox       = [1 1 1]; % here is the voxel size
    flags.interp    = 0; % 0: nearest neighbor, 4: 4th degree b spline
    flags.wrap      = [0 0 0];
    flags.prefix    = 'w';
    
    % normalize the image with PC1 weights
    job.subj.vol{1} = niAnatomy;
    job.subj.resample{1} = pc1_NameSave;
    job.woptions = flags;
    spm_run_norm(job);
    
    % normalize the image with PC2 weights
    job.subj.vol{1} = niAnatomy;
    job.subj.resample{1} = pc2_NameSave;
    job.woptions = flags;
    spm_run_norm(job);
    
    % normalize the image with slope weights
    job.subj.vol{1} = niAnatomy;
    job.subj.resample{1} = slope_NameSave;
    job.woptions = flags;
    spm_run_norm(job);
    
    job.subj.vol{1} = niAnatomy;
    job.subj.resample{1} = onset_NameSave;
    job.woptions = flags;
    spm_run_norm(job);
    
    job.subj.vol{1} = niAnatomy;
    job.subj.resample{1} = fwhm_NameSave;
    job.woptions = flags;
    spm_run_norm(job);
end

%% load all subjects and add in one MNI image 

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,1,1,1,1};

flags.bb    = [-90 -120 -60; 90 96 130];
all_mni_pc1 = NaN([diff(flags.bb)+1 length(sub_labels)]);
all_mni_pc2 = NaN([diff(flags.bb)+1 length(sub_labels)]);
all_mni_cod = NaN([diff(flags.bb)+1 length(sub_labels)]);
all_mni_slope = NaN([diff(flags.bb)+1 length(sub_labels)]);
all_mni_onset = NaN([diff(flags.bb)+1 length(sub_labels)]);
all_mni_fwhm = NaN([diff(flags.bb)+1 length(sub_labels)]);

data_in = 'PPG';
for ss = 1:length(sub_labels)

    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);

    % get PC1, PC2, slope, onset, fwhm
    pc1_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wsub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_canoPc1Weights.nii']));
    pc2_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wsub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_canoPc2Weights.nii']));
    slope_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wsub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_modelSlope.nii']));
    onset_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wsub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_modelOnset.nii']));
    fwhm_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wsub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_model_FWHM.nii']));
    
    all_mni_pc1(:,:,:,ss) = pc1_mni.data;
    all_mni_pc2(:,:,:,ss) = pc2_mni.data;
    all_mni_slope(:,:,:,ss) = slope_mni.data;
    all_mni_onset(:,:,:,ss) = onset_mni.data;
    all_mni_fwhm(:,:,:,ss) = fwhm_mni.data;
    
    % get COD to set a reliability threshold
    wfcod = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wfsub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_codPPG.nii']));
        
    all_mni_cod(:,:,:,ss) = wfcod.data;
end

% Slice view of beta1 (pc1) and beta2 (pc2) using fancy color circle
niAnatomy = niftiRead(fullfile(dDir,'derivatives','mni','rsingle_subj_T1.nii'));
imDims = [-90 -120 -70; 90 100 90];
curPos = [1 12 4]; 
acpcXform = pc1_mni.qto_xyz;

%% Figure 6A panels: timing of local minima with onset for color, R^2 for intensity
onset_mean = mean(all_mni_onset,4); 
slope_mean = mean(all_mni_slope,4); 
cod_mean = mean(all_mni_cod,4); 

% plot blood
onset_mean(slope_mean>0) = NaN;
cod_mean(slope_mean>0) = NaN;

onset_plot = pc1_mni; % initialize nifti
cod_plot = pc1_mni; % initialize nifti
onset_plot.data = onset_mean;
cod_plot.data = cod_mean;

% plot Saggital figure
sliceThisDim = 1;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,2);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['Figure6A_MNI_Saggital' int2str(curPos(sliceThisDim)) '_vneg2']))

% Plot Axial
sliceThisDim = 3;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,2);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['Figure6A_MNI_Axial' int2str(curPos(sliceThisDim)) '_vneg2']))

% Plot Coronal
sliceThisDim = 2;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,2);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['Figure6A_MNI_Coronal' int2str(curPos(sliceThisDim)) '_vneg2']))

%% Figure 6B panels: timing of local maxima with onset for color, R^2 for intensity

onset_mean = mean(all_mni_onset,4); 
slope_mean = mean(all_mni_slope,4); 
cod_mean = mean(all_mni_cod,4); 

% plot csf
onset_mean(slope_mean<0) = NaN;
cod_mean(slope_mean<0) = NaN;

onset_plot = pc1_mni;
cod_plot = pc1_mni;
onset_plot.data = onset_mean;
cod_plot.data = cod_mean;
acpcXform = pc1_mni.qto_xyz;

% plot Saggital figure
sliceThisDim = 1;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,3);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['Figure6B_MNI_Saggital' int2str(curPos(sliceThisDim)) '_vpos2']))

% Plot Axial
sliceThisDim = 3;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,3);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['Figure6B_MNI_Axial' int2str(curPos(sliceThisDim)) '_vpos2']))

% Plot Coronal
sliceThisDim = 2;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,3);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['Figure6B_MNI_Coronal' int2str(curPos(sliceThisDim)) '_vpos2']))


%%
%% Uncomment below if you want to make figures of the width distribution
%%
% 
% % plot width of local minima 
% width_mean = mean(all_mni_fwhm,4); 
% slope_mean = mean(all_mni_slope,4); 
% cod_mean = mean(all_mni_cod,4); 
% 
% % plot blood
% width_mean(slope_mean>0) = NaN;
% cod_mean(slope_mean>0) = NaN;
% 
% width_plot = pc1_mni; % initialize nifti
% cod_plot = pc1_mni; % initialize nifti
% width_plot.data = width_mean;
% cod_plot.data = cod_mean;
% 
% % plot Saggital figure
% sliceThisDim = 1;
% bbOverlayDotsAnat_Color2D(width_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,1);
% set(gca,'CLim',[0 1])
% 
% 
% % plot width of local maxima
% width_mean = mean(all_mni_fwhm,4); 
% slope_mean = mean(all_mni_slope,4); 
% cod_mean = mean(all_mni_cod,4); 
% 
% % plot CSF
% width_mean(slope_mean<0) = NaN;
% cod_mean(slope_mean<0) = NaN;
% 
% width_plot = pc1_mni; % initialize nifti
% cod_plot = pc1_mni; % initialize nifti
% width_plot.data = width_mean;
% cod_plot.data = cod_mean;
% 
% % plot Saggital figure
% sliceThisDim = 1;
% bbOverlayDotsAnat_Color2D(width_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,1);
% set(gca,'CLim',[0 1])


