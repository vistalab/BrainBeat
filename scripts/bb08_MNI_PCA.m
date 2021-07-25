

%% The functionals were not writted in the same space as the T1 
% here save the coregistration matrix to the f.._codPPG.nii in sto_xyz/ijk
% this is necessary for the normalization

% dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

%%%% TO CHANGE: pc1 and pc2 are now saved in T1-space, update below

data_in = 'PPG';
s = 6;
scan_nr = 2;

s_info = bb_subs(s);
subj = s_info.subj;
    
scan = s_info.scan{scan_nr};
scanName = s_info.scanName{scan_nr};

% name for pc1 and pc2 beta weights
pc1_Name = fullfile(dDir,subj,scan,[scanName '_' data_in '_pc1.nii.gz']);
pc2_Name = fullfile(dDir,subj,scan,[scanName '_' data_in '_pc2.nii.gz']);
% read
pc1 = niftiRead(pc1_Name);
pc2 = niftiRead(pc2_Name);

% load coregistration matrix for the functionals
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% save this coregistration matrix in the f_...codPPG.nii
pc1.qto_xyz = acpcXform;
pc1.qto_ijk = inv(acpcXform);
pc1.sto_xyz = acpcXform;
pc1.sto_ijk = inv(acpcXform);

pc2.qto_xyz = acpcXform;
pc2.qto_ijk = inv(acpcXform);
pc2.sto_xyz = acpcXform;
pc2.sto_ijk = inv(acpcXform);

% name for pc1 & pc2
pc1_NameSave = fullfile(dDir,subj,scan,['f' scanName '_' data_in '_pc1w.nii']);
pc2_NameSave = fullfile(dDir,subj,scan,['f' scanName '_' data_in '_pc2w.nii']);

% save pca1 & pc2
niftiWrite(pc1,pc1_NameSave);
niftiWrite(pc2,pc2_NameSave);


%%%%%% now we can normalize the PC1 and PC2 beta weight image

spm('Defaults','fmri')

% A) volume where parameters are estimated: anat
% code will search for y_anat.nii file to use for normalization
niAnatomy = fullfile(dDir,subj,s_info.anat,['f' s_info.anatName '.nii']);

% B) volume you want to normalize, needs to be coregistered with A)
% code will append a w before this file
pc1_NameSave = fullfile(dDir,subj,scan,['f' scanName '_' data_in '_pc1w.nii']);
pc2_NameSave = fullfile(dDir,subj,scan,['f' scanName '_' data_in '_pc2w.nii']);

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

subj_inds   = [1 2 3 4 5 6];
scan_inds   = [3 3 3 1 1 2];
flags.bb    = [-90 -120 -60; 90 96 130];
all_mni_pc1 = NaN([diff(flags.bb)+1 length(subj_inds)]);
all_mni_pc2 = NaN([diff(flags.bb)+1 length(subj_inds)]);
all_mni_cod = NaN([diff(flags.bb)+1 length(subj_inds)]);

data_in = 'PPG';
for aa = 1:length(subj_inds)

    s = subj_inds(aa);
    scan_nr = scan_inds(aa);

    s_info = bb_subs(s);
    subj = s_info.subj;

    scan = s_info.scan{scan_nr};
    scanName = s_info.scanName{scan_nr};

    pc1_mni = niftiRead(fullfile(dDir,subj,scan,['wf' scanName '_' data_in '_pc1w.nii']));
    pc2_mni = niftiRead(fullfile(dDir,subj,scan,['wf' scanName '_' data_in '_pc2w.nii']));
    
    all_mni_pc1(:,:,:,aa) = pc1_mni.data;
    all_mni_pc2(:,:,:,aa) = pc2_mni.data;
    
    % also get COD to set a reliability threshold
    wfcod = niftiRead(fullfile(dDir,subj,scan, ['wf' scanName '_codPPG.nii']));

    all_mni_cod(:,:,:,aa) = wfcod.data;
end

%% 

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
% working on this

Nthreshold = 4;
pc1_posneg = {'pc1pos','pc1neg'};

for plot_positive = 1:2 % 1 = positive, 2 = negative, 

    if plot_positive==1
        select_voxels = find(pc1_mni.data>=Nthreshold);
    elseif plot_positive==2
        select_voxels = find(pc1_mni.data<=-Nthreshold);
    end

    % Get indiced of selected voxels
    [ii,jj,kk] = ind2sub(size(pc1_mni.data),select_voxels);
    ijk_func = [ii jj kk];
    clear ii jj kk % housekeeping

    % average PC 2 for color
    pc2_mean = mean(all_mni_pc2,4);
    ijk_color = pc2_mean(select_voxels);
    ijk_colorInd = (ijk_color-min(ijk_color));
    ijk_colorInd = 1+round(99*ijk_colorInd/max(ijk_colorInd));

    if plot_positive==1
        cm = jet(100);
        cm = cm(end:-1:1,:);
    elseif plot_positive==2
        cm = jet(100);
    end

    % Get mni coordinates of voxels 
    xyz_mni = mrAnatXformCoords(wfcod.sto_xyz, ijk_func);


    % load MNI rendings 
    load(fullfile(dDir,'MNI_cortex_left.mat'))
    gl.vertices = cortex.vert;
    gl.faces = cortex.tri;
    gl.mat = [1 0 0 1;0 1 0 1; 0 0 1 1; 0 0 0 1];
    gl = gifti(gl);

    load(fullfile(dDir,'MNI_cortex_right.mat'))
    gr.vertices = cortex.vert;
    gr.faces = cortex.tri;
    gr.mat = [1 0 0 1;0 1 0 1; 0 0 1 1; 0 0 0 1];
    gr = gifti(gr);

    % plot right
    figure,hold on
    ieeg_RenderGifti(gr)
    % add PC2 in color:
    for kk = 1:size(xyz_mni,1)
        if xyz_mni(kk,1)>-10
            plot3(xyz_mni(kk,1),xyz_mni(kk,2),xyz_mni(kk,3),'.','Color',cm(ijk_colorInd(kk),:))
        end
    end
    set(gcf,'PaperPositionMode','auto') 
    ieeg_viewLight(270,0)
    print('-painters','-r300','-dpng',[dDir '/figures/reliable/mni_all_right_render' pc1_posneg{plot_positive} '_view1_v00'])
    ieeg_viewLight(90,0)
    print('-painters','-r300','-dpng',[dDir '/figures/reliable/mni_all_right_render' pc1_posneg{plot_positive} '_view2_v00'])
    
    % plot left
    figure,hold on
    ieeg_RenderGifti(gl)
    % add PC2 in color:
    for kk = 1:size(xyz_mni,1)
        if xyz_mni(kk,1)<10
            plot3(xyz_mni(kk,1),xyz_mni(kk,2),xyz_mni(kk,3),'.','Color',cm(ijk_colorInd(kk),:))
        end
    end
    set(gcf,'PaperPositionMode','auto') 
    ieeg_viewLight(270,0)
    print('-painters','-r300','-dpng',[dDir '/figures/reliable/mni_all_left_render' pc1_posneg{plot_positive} '_view1_v00'])
    ieeg_viewLight(90,0)
    print('-painters','-r300','-dpng',[dDir '/figures/reliable/mni_all_left_render' pc1_posneg{plot_positive} '_view2_v00'])

end

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



