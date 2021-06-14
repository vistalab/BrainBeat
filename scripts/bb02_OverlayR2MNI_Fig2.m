
% make sure full path is added
spm('Defaults','fmri')
spm_path = fileparts(which('spm'));
addpath(fullfile(spm_path,'config'))

%% The functionals were not writted in the same space as the T1 
% here save the coregistration matrix to the f.._codPPG.nii in sto_xyz/ijk
% this is necessary for the normalization

% dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

sub_labels = {'5'}; 
ses_labels = {'1'}; 
acq_labels = {'4mmFA48'};
run_nrs = {[1]};

ss = 1;%:length(sub_labels) % subjects/ses/acq
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
save_name_base = (['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

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

%% load all subjects and render MNI

Rthreshold = 0.9;
flags.bb    = [-90 -120 -60; 90 96 130];

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],1,1,1,1,1};

% initialize output for COD in MNI space for all subjects
all_mni_cod = NaN([diff(flags.bb)+1 length(sub_labels)]);

% just take one run
rr = 1;% run_nr

all_mni_cod = [];

for ss = 1:length(sub_labels)
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};

    run_nr = run_nrs{ss}(rr);

    save_dir = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label]);
    save_name_base = (['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

    wfcod = niftiRead(fullfile(save_dir,['wf' save_name_base '_codPPG.nii']));

    % get COD in MNI space for all subjects
    all_mni_cod(:,:,:,ss) = wfcod.data;
    
    % Get mni coordinates of voxels above threshold
    select_voxels = find(wfcod.data>=Rthreshold);

    % Get indiced of selected voxels
    [ii,jj,kk] = ind2sub(size(wfcod.data),select_voxels);
    ijk_func = [ii jj kk];
    clear ii jj kk % housekeeping

    all_mni(ss).xyz_mni = mrAnatXformCoords(wfcod.sto_xyz, ijk_func);
end

%% 

save_dir = fullfile(dDir,'derivatives','brainbeat','group');
save_name = fullfile(save_dir,['group_4mmFA48_wfcodPPG.nii']);

cod_mni = wfcod;

cod_mni.data = sum(all_mni_cod>0.7,4);

niftiWrite(cod_mni,save_name);

%%
%% left off here
%%
figure,hold on

for aa = 1:length(subj_inds)
    x = all_mni(aa).xyz_mni(:,1);
    y = all_mni(aa).xyz_mni(:,2);
    z = all_mni(aa).xyz_mni(:,3);
    plot3(x,y,z,'b.')
end

%% test with rendering

load(fullfile(dDir,'MNI_cortex_left.mat'))

g.vertices = cortex.vert;
g.faces = cortex.tri;
g.mat = [1 0 0 1;0 1 0 1; 0 0 1 1; 0 0 0 1];
g = gifti(g);

p = [all_mni(1).xyz_mni; all_mni(2).xyz_mni; all_mni(3).xyz_mni; all_mni(4).xyz_mni; all_mni(5).xyz_mni; all_mni(6).xyz_mni];

this_p = p;
this_p(p(:,1)>30,:) = [];
figure
ieeg_RenderGifti(g)
ieeg_elAdd(this_p,[.5 0 1],5)

% for aa = 1:length(subj_inds)
%     select_these_point = all_mni(aa).xyz_mniloc_view(:,1)<10; % left
%     x = all_mni(aa).xyz_mni(select_these_point,1);
%     y = all_mni(aa).xyz_mni(select_these_point,2);
%     z = all_mni(aa).xyz_mni(select_these_point,3);
%     scatter3(x,y,z,'filled','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
% end

set(gcf,'PaperPositionMode','auto') 
ieeg_viewLight(270,0)
print('-painters','-r300','-dpng',[dDir '/figures/reliable/mni_left_render_V1_codth0_9'])
ieeg_viewLight(90,0)
print('-painters','-r300','-dpng',[dDir '/figures/reliable/mni_left_render_V2_codth0_9'])

    
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



