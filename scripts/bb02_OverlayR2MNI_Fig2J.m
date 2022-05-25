

% make sure full path is added
spm('Defaults','fmri')
spm_path = fileparts(which('spm'));
addpath(fullfile(spm_path,'config'))

[~,dDir] = bbPath;

%% The functionals were not writted in the same space as the T1 
% here save the coregistration matrix to the f.._codPPG.nii in sto_xyz/ijk
% this is necessary for the normalization

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,1,1,1,1};

for ss = 1:length(sub_labels) % subjects/ses/acq
    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    
    run_nr = run_nrs{ss}(rr);
    
    % Get the anatomy used for SPM segmentation:
    t1w_BIDSnameSPM = fullfile('derivatives','spmSegmentation',['sub-' sub_label],['ses-' ses_label],...
                ['sub-' sub_label '_ses-' ses_label '_T1w.nii']);
    niAnatomySPM = niftiRead(fullfile(dDir,t1w_BIDSnameSPM));
    
    save_dir = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label]);
    save_name_base = (['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);
    
    % load coregistration matrix (for the functionals):
    load(fullfile(save_dir,[save_name_base '_AcpcXform_new.mat']));
    acpcXform = acpcXform_new; clear acpcXform_new
    
    % reliability (COD)
    ppgRname = fullfile(save_dir,[save_name_base '_codPPG.nii.gz']);
    ppgR = niftiRead(ppgRname); % correlation with PPG
    ppgRnameT1spaceName = fullfile(save_dir,['f' save_name_base '_codPPG.nii']);
    
    % save COD in T1w space
    ni  = ppgR;
    % save this coregistration matrix in the f_...codPPG.nii 
    ni.qto_xyz = acpcXform;
    ni.qto_ijk = inv(acpcXform);
    ni.sto_xyz = acpcXform;
    ni.sto_ijk = inv(acpcXform);
    
    niftiWrite(ni,ppgRnameT1spaceName);
    
    %%%%%% now we normalize the new fcod image
    
    % A) volume where parameters are estimated: anat
    % code will search for y_anat.nii file to use for normalization
    % t1w_BIDSnameSPM
    
    % B) volume you want to normalize, needs to be coregistered with A)
    % code will append a w before this file
    % ppgRnameT1spaceName;
    
    flags.preserve  = 0;
    flags.bb        = [-90 -120 -60; 90 96 130];
    flags.vox       = [1 1 1]; % here is the voxel size
    flags.interp    = 0; % 0: neirest neighbor, 4: 4th degree b spline
    flags.wrap      = [0 0 0];
    flags.prefix    = 'w';
    
    job.subj.vol{1} = fullfile(dDir,t1w_BIDSnameSPM);
    job.subj.resample{1} = ppgRnameT1spaceName;
    job.woptions = flags;
    
    % normalize the image with COD
    spm_run_norm(job);
end

%%
%% load all subjects and add in one MNI image 

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

flags.bb    = [-90 -120 -60; 90 96 130];
all_mni_cod = NaN([diff(flags.bb)+1 length(sub_labels)]);

data_in = 'PPG';
for ss = 1:length(sub_labels)

    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);
    
    % get COD to set a reliability threshold
    wfcod = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wfsub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_codPPG.nii']));
        
    all_mni_cod(:,:,:,ss) = wfcod.data;
end


%%
%% Slice view group R
%%

% Slice view of beta1 (pc1) and beta2 (pc2) using fancy color circle
niAnatomy = niftiRead(fullfile(dDir,'derivatives','mni','rsingle_subj_T1.nii'));
imDims = [-90 -120 -70; 90 100 90];
curPos = [1 12 4]; 

% Onset for color, cod for intensity
cod_mean = mean(all_mni_cod,4); 
cod_mean(cod_mean<.3) = 0; 

load loc_colormap
% cm = parula(100);

% get structure and xform for a singel subjects 
cod_plot = wfcod; % load single subject just for structure.
cod_plot.data = cod_mean; % replace with group average
acpcXform = wfcod.qto_xyz; 

set(gca,'CLim',[0 1])

% plot Saggital figure
sliceThisDim = 1;
figure
bbOverlayDotsAnat_PickColor(cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,cm(33:end,:),1)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['Figure2J_MNI_Saggital' int2str(curPos(sliceThisDim)) '_vCOD']))

% Plot Axial
sliceThisDim = 3;
figure
bbOverlayDotsAnat_PickColor(cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,cm(33:end,:),1)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['Figure2J_MNI_Axial' int2str(curPos(sliceThisDim)) '_vCOD']))

% Plot Coronal
sliceThisDim = 2;
figure
bbOverlayDotsAnat_PickColor(cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,cm(33:end,:),1)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['Figure2J_MNI_Coronal' int2str(curPos(sliceThisDim)) '_vCOD']))







