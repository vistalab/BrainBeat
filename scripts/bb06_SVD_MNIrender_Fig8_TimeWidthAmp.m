

%% Normalize PC1/PC2 and the slope, latency and width

% dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

ss = 6;
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
% niAnatomy = niftiRead(fullfile(dDir,t1w_BIDSname));

% B) volume you want to normalize, needs to be coregistered with A)
% code will append a w before this file
pc1_NameSave = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['fsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_canoPc1Weights.nii']);
pc2_NameSave = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['fsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_canoPc2Weights.nii']);

% Get base name
slope_NameSave = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['fsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_modelSlope.nii']);
onset_NameSave = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['fsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_modelOnset.nii']);
fwhm_NameSave = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['fsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_model_FWHM.nii']);

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

%% load all subjects and add in one MNI image 

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

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
        ['wfsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_canoPc1Weights.nii']));
    pc2_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wfsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_canoPc2Weights.nii']));
    slope_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wfsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_modelSlope.nii']));
    onset_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wfsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_modelOnset.nii']));
    fwhm_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wfsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_model_FWHM.nii']));
    
    all_mni_pc1(:,:,:,ss) = pc1_mni.data;
    all_mni_pc2(:,:,:,ss) = pc2_mni.data;
    all_mni_slope(:,:,:,ss) = slope_mni.data;
    all_mni_onset(:,:,:,ss) = onset_mni.data;
    all_mni_fwhm(:,:,:,ss) = fwhm_mni.data;
    
    % get COD to set a reliability threshold
    wfcod = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wfsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_codPPG.nii']));
        
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

cod_plot = pc1_mni;
cod_plot.data = cod_mean;

set(gca,'CLim',[0 1])

% plot Saggital figure
sliceThisDim = 1;
figure
bbOverlayDotsAnat_PickColor(cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,cm(33:end,:),1)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['MNI_Saggital' int2str(curPos(sliceThisDim)) '_vCOD']))

% Plot Axial
sliceThisDim = 3;
figure
bbOverlayDotsAnat_PickColor(cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,cm(33:end,:),1)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['MNI_Axial' int2str(curPos(sliceThisDim)) '_vCOD']))

% Plot Coronal
sliceThisDim = 2;
figure
bbOverlayDotsAnat_PickColor(cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,cm(33:end,:),1)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['MNI_Coronal' int2str(curPos(sliceThisDim)) '_vCOD']))

%% Negative peaks: onset for color, cod for intensity
onset_mean = mean(all_mni_onset,4); 
slope_mean = mean(all_mni_slope,4); 
cod_mean = mean(all_mni_cod,4); 

% plot blood
onset_mean(slope_mean>0) = NaN;
cod_mean(slope_mean>0) = NaN;

onset_plot = pc1_mni;
cod_plot = pc1_mni;
onset_plot.data = onset_mean;
cod_plot.data = cod_mean;
acpcXform = pc1_mni.qto_xyz;

curPos = [1 12 4]; 

% plot Saggital figure
sliceThisDim = 1;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,2);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['MNI_Saggital' int2str(curPos(sliceThisDim)) '_vneg2']))

% Plot Axial
sliceThisDim = 3;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,2);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['MNI_Axial' int2str(curPos(sliceThisDim)) '_vneg2']))

% Plot Coronal
sliceThisDim = 2;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,2);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['MNI_Coronal' int2str(curPos(sliceThisDim)) '_vneg2']))

%% Positive peaks: onset for color, cod for intensity
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

curPos = [1 12 4]; 

% plot Saggital figure
sliceThisDim = 1;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,3);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['MNI_Saggital' int2str(curPos(sliceThisDim)) '_vpos2']))

% Plot Axial
sliceThisDim = 3;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,3);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['MNI_Axial' int2str(curPos(sliceThisDim)) '_vpos2']))

% Plot Coronal
sliceThisDim = 2;
bbOverlayDotsAnat_Color2D(onset_plot,cod_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .3],8,.3,3);
set(gca,'CLim',[0 1])
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['MNI_Coronal' int2str(curPos(sliceThisDim)) '_vpos2']))



%% 
%%
%% Renderings code not implemented yet - takes take a few minutes for each
%%
%%

Rthreshold = 0.5;
these_areas = all_mni_cod>=Rthreshold;
these_areas = sum(these_areas,4);
select_voxels = find(these_areas>=3); % voxel>Rthreshold in more than 3 subjects

% Get indiced of selected voxels
[ii,jj,kk] = ind2sub(size(wfcod.data),select_voxels);
ijk_func = [ii jj kk];
clear ii jj kk % housekeeping

% Get mni coordinates of voxels
xyz_mni = mrAnatXformCoords(wfcod.sto_xyz, ijk_func);


%% % load MNI cortical surfaces 
load(fullfile(dDir,'derivatives','mni','MNI_cortex_left.mat'))
gl.vertices = cortex.vert;
gl.faces = cortex.tri;
gl.mat = [1 0 0 1;0 1 0 1; 0 0 1 1; 0 0 0 1];
gl = gifti(gl);

load(fullfile(dDir,'derivatives','mni','MNI_cortex_right.mat'))
gr.vertices = cortex.vert;
gr.faces = cortex.tri;
gr.mat = [1 0 0 1;0 1 0 1; 0 0 1 1; 0 0 0 1];
gr = gifti(gr);

%% plot right
figure,hold on
fg = ieeg_RenderGifti(gr);
% add PC2 in color:
for kk = 1:size(xyz_mni,1)
    if xyz_mni(kk,1)>-10
        plot3(xyz_mni(kk,1),xyz_mni(kk,2),xyz_mni(kk,3),'.','Color',data_colors_rgb(kk,:))
    end
end
set(fg,'FaceAlpha',.5)
set(gcf,'PaperPositionMode','auto') 
ieeg_viewLight(270,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_right_render_view1_v02'))

ieeg_viewLight(90,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_right_render_view2_v02'))

%% plot left
figure,hold on
fg = ieeg_RenderGifti(gl);
% add PC2 in color:
for kk = 1:size(xyz_mni,1)
    if xyz_mni(kk,1)<10
        plot3(xyz_mni(kk,1),xyz_mni(kk,2),xyz_mni(kk,3),'.','Color',data_colors_rgb(kk,:))
    end
end
set(fg,'FaceAlpha',.5)
set(gcf,'PaperPositionMode','auto') 
ieeg_viewLight(270,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_left_render_view1_v02'))
ieeg_viewLight(90,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_left_render_view2_v02'))



%%
%% test only plotting positive/negative
%%

pc_complex = complex(pc1_mean(select_voxels),pc2_mean(select_voxels));
pc_angle = angle(pc_complex*(1-1i)); % multiply by (1-i) to rotatio 45 deg

%% plot right
figure,hold on
fg = ieeg_RenderGifti(gr);
% add PC2 in color:
for kk = 1:size(xyz_mni,1)
    if xyz_mni(kk,1)>-10
        if pc_angle(kk)<0
            plot3(xyz_mni(kk,1),xyz_mni(kk,2),xyz_mni(kk,3),'.','Color',data_colors_rgb(kk,:))
        end
    end
end
set(fg,'FaceAlpha',.5)
set(gcf,'PaperPositionMode','auto') 
ieeg_viewLight(270,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_right_render_view1_vPos'))

ieeg_viewLight(90,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_right_render_view2_vPos'))

%% plot left
figure,hold on
fg = ieeg_RenderGifti(gl);
% add PC2 in color:
for kk = 1:size(xyz_mni,1)
    if xyz_mni(kk,1)<10
        if pc_angle(kk)<0
            plot3(xyz_mni(kk,1),xyz_mni(kk,2),xyz_mni(kk,3),'.','Color',data_colors_rgb(kk,:))
        end
    end
end

set(fg,'FaceAlpha',.5)
set(gcf,'PaperPositionMode','auto') 
ieeg_viewLight(270,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_left_render_view1_vPos'))
ieeg_viewLight(90,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_left_render_view2_vPos'))
