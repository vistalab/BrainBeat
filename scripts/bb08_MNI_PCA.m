

%% The functionals were not writted in the same space as the T1 
% here save the coregistration matrix to the f.._codPPG.nii in sto_xyz/ijk
% this is necessary for the normalization

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


%% load all subjects and add in one MNI image 

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

flags.bb    = [-90 -120 -60; 90 96 130];
all_mni_pc1 = NaN([diff(flags.bb)+1 length(sub_labels)]);
all_mni_pc2 = NaN([diff(flags.bb)+1 length(sub_labels)]);
all_mni_cod = NaN([diff(flags.bb)+1 length(sub_labels)]);

data_in = 'PPG';
for ss = 1:length(sub_labels)

    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);

    % get PC1 and PC2 
    pc1_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wfsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_canoPc1Weights.nii']));
    pc2_mni = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wfsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_space-T1w_canoPc2Weights.nii']));
    
    all_mni_pc1(:,:,:,ss) = pc1_mni.data;
    all_mni_pc2(:,:,:,ss) = pc2_mni.data;
    
    % get COD to set a reliability threshold
    wfcod = niftiRead(fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['wfsub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_codPPG.nii']));
        
    all_mni_cod(:,:,:,ss) = wfcod.data;
end

%% did not do this step yet

% name for pc1 & pc2
pc1_MNI_all = './local/allMNI_pc1w.nii.gz';
pc2_MNI_all = './local/allMNI_pc2w.nii.gz';
cod_MNI_all = './local/allMNI_cod.nii.gz';

% save pc1 beta weights, threshold by cod, and add across subjects
pc1_mni.data = all_mni_pc1;
cod_th = 0.1;
pc1_mni.data(all_mni_cod<cod_th) = 0;
pc1_mni.data(all_mni_cod>=cod_th & all_mni_pc1>0) = 1; % 1 to everything pc1>0
pc1_mni.data(all_mni_cod>=cod_th & all_mni_pc1<0) = -1; % -1 to everything pc1<0
pc1_mni.data = sum(pc1_mni.data,4); % sum across subjects
niftiWrite(pc1_mni,pc1_MNI_all);


%% LEFT OFF HERE, not writing PC2 yet
% pc2_mni.data = all_mni_pc2;
% niftiWrite(pc2_mni,pc2_MNI_all);
% 
% cod_mni = pc1_mni;
% cod_mni.data = all_mni_cod;
% niftiWrite(cod_mni,cod_MNI_all);


%%
%% load all subjects and render MNI
%%

Rthreshold = 0.7;
select_voxels = find(median(all_mni_cod(:,:,:,1:5),4)>=Rthreshold);

% Get indiced of selected voxels
[ii,jj,kk] = ind2sub(size(wfcod.data),select_voxels);
ijk_func = [ii jj kk];
clear ii jj kk % housekeeping

% average PC1 and PC2 for color
pc1_mean = mean(all_mni_pc1,4);
pc2_mean = mean(all_mni_pc2,4);
pc1_mean(pc1_mean<-1) = -1;
pc2_mean(pc2_mean<-1) = -1;
pc1_mean(pc1_mean>1) = 1;
pc2_mean(pc2_mean>1) = 1;
data_colors_rgb = bbData2Colors([pc1_mean(select_voxels) pc2_mean(select_voxels)]);
    
% Get mni coordinates of voxels 
xyz_mni = mrAnatXformCoords(wfcod.sto_xyz, ijk_func);

% load MNI rendings 
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

% plot right
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
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_right_render_view1_v01'))
ieeg_viewLight(90,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_right_render_view2_v01'))

% plot left
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
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_left_render_view1_v01'))
ieeg_viewLight(90,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_left_render_view2_v01'))

%%
%% now try to render/combine all subjects above trheshold
%%

Rthreshold = 0.8;

% load MNI rendings 
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


% plot right
figure,hold on
fg = ieeg_RenderGifti(gr);

for ss = 1:5
    select_voxels = find(all_mni_cod(:,:,:,ss)>=Rthreshold);

    % Get indiced of selected voxels
    [ii,jj,kk] = ind2sub(size(wfcod.data),select_voxels);
    ijk_func = [ii jj kk];
    clear ii jj kk % housekeeping

    % average PC1 and PC2 for color
    pc1_mean = all_mni_pc1(:,:,:,ss);
    pc2_mean = all_mni_pc2(:,:,:,ss);
    pc1_mean(pc1_mean<-1) = -1;
    pc2_mean(pc2_mean<-1) = -1;
    pc1_mean(pc1_mean>1) = 1;
    pc2_mean(pc2_mean>1) = 1;
    data_colors_rgb = bbData2Colors([pc1_mean(select_voxels) pc2_mean(select_voxels)]);

    % Get mni coordinates of voxels 
    xyz_mni = mrAnatXformCoords(wfcod.sto_xyz, ijk_func);

    % add PC2 in color:
    for kk = 1:size(xyz_mni,1)
        if xyz_mni(kk,1)>-10
            plot3(xyz_mni(kk,1),xyz_mni(kk,2),xyz_mni(kk,3),'.','Color',data_colors_rgb(kk,:))
        end
    end
end
set(fg,'FaceAlpha',.5)
set(gcf,'PaperPositionMode','auto') 
ieeg_viewLight(270,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_right_render_view1_v02'))
ieeg_viewLight(90,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_right_render_view2_v02'))


% plot left
figure,hold on
fg = ieeg_RenderGifti(gl);

for ss = 1:5
    select_voxels = find(all_mni_cod(:,:,:,ss)>=Rthreshold);

    % Get indiced of selected voxels
    [ii,jj,kk] = ind2sub(size(wfcod.data),select_voxels);
    ijk_func = [ii jj kk];
    clear ii jj kk % housekeeping

    % average PC1 and PC2 for color
    pc1_mean = all_mni_pc1(:,:,:,ss);
    pc2_mean = all_mni_pc2(:,:,:,ss);
    pc1_mean(pc1_mean<-1) = -1;
    pc2_mean(pc2_mean<-1) = -1;
    pc1_mean(pc1_mean>1) = 1;
    pc2_mean(pc2_mean>1) = 1;
    data_colors_rgb = bbData2Colors([pc1_mean(select_voxels) pc2_mean(select_voxels)]);

    % Get mni coordinates of voxels 
    xyz_mni = mrAnatXformCoords(wfcod.sto_xyz, ijk_func);

    for kk = 1:size(xyz_mni,1)
        if xyz_mni(kk,1)<10
            plot3(xyz_mni(kk,1),xyz_mni(kk,2),xyz_mni(kk,3),'.','Color',data_colors_rgb(kk,:))
        end
        
    end
end
set(fg,'FaceAlpha',.5)
set(gcf,'PaperPositionMode','auto') 
ieeg_viewLight(270,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_left_render_view1_v02'))
ieeg_viewLight(90,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','mni_all_left_render_view2_v02'))

close all

%%
%%
%%

%% test
% addpath to MyCrust270110
% p = [all_mni(1).xyz_mni; all_mni(2).xyz_mni; all_mni(3).xyz_mni; all_mni(4).xyz_mni; all_mni(5).xyz_mni; all_mni(6).xyz_mni];
p = [all_mni(1).xyz_mni; all_mni(2).xyz_mni; all_mni(3).xyz_mni];


%% Run  program
[t]=MyRobustCrust(p);


%% plot the points cloud
figure(1);
set(gcf,'position',[0,0,1280,800]);
subplot(1,2,1)
hold on
axis equal
title('Points Cloud','fontsize',14)
plot3(p(:,1),p(:,2),p(:,3),'g.')
axis vis3d
view(3)


%% plot the output triangulation
figure(1)
subplot(1,2,2)
hold on
title('Output Triangulation','fontsize',14)
axis equal
trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')%plot della superficie
axis vis3d
view(3)



