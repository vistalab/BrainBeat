
% This script generates Figure 6 from the manuscript titled:
%
%%%%% Measuring brain beats: cardiac-aligned fast fMRI signals %%%%%
% Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell
% 

%%
%% FIGURES & RENDERINGS in individual subjects
%%

clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%% Plot R/Betas back on brain

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

ss = 1;
rr = 1;% run_nr
sub_label = sub_labels{ss};
ses_label = ses_labels{ss};
acq_label = acq_labels{ss};
run_nr = run_nrs{ss}(rr);

%% Get base name of saved data

save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);

% load pc1
pc1Weight = niftiRead([save_name_base '_space-T1w_canoPc1Weights.nii.gz']);
% load pc2
pc2Weight = niftiRead([save_name_base '_space-T1w_canoPc2Weights.nii.gz']);
% load r weights (how much explained by canonical PC1-2)
r_weight = niftiRead([save_name_base '_space-T1w_canoPc12R.nii.gz']);
acpcXform = pc1Weight.qto_xyz;

% Get anatomy
t1w_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'anat',...
            ['sub-' sub_label '_ses-' ses_label '_T1w.nii.gz']);
niAnatomy = niftiRead(fullfile(dDir,t1w_BIDSname));


% plot full circle as in supplementary figure S4
if ss == 1
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-6 26 17]; 
elseif ss == 2
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [-1 50 -21]; 
elseif ss == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [-2 26 -63]; 
elseif ss == 4
    imDims = [-90 -120 -50; 90 130 120];
    curPos = [0 4 35];
elseif ss == 5
    imDims = [-90 -120 -100; 90 130 120];
    curPos = [-2,18,38]; 
elseif ss == 6
    imDims = [-90 -120 -100; 90 130 120];
    curPos = [-4,18,38];
end

% plot beta1 (pc1) and beta2 (pc2) using fancy color circle
sliceThisDim = 1;
bbOverlayDotsAnat_FancyColorCircle(pc1Weight,pc2Weight,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.7);
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['subj' int2str(ss) '_run' int2str(rr) '_exampleAxial38']))

% Get functional voxels to plot with rendering: brain mask and COD
% threshold
%%%% get a segmentation to create a mask for the brain, CSF and vessels:
% niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
niSegm2 = niftiRead([save_name_base '_combineSegm.nii.gz']);
brainMask = niSegm2.data>0;
%%%% load PPG-RCOD
Rthreshold = .5; % only plot reliable voxels
ppgR = niftiRead([save_name_base '_codPPG.nii.gz']);
%%%% Use mask and R threshold:
select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);

% PC1/PC2 weights for voxels to plot:
pc12_render = [pc1Weight.data(select_voxels) pc2Weight.data(select_voxels)];

% Get indiced of selected voxels
[ii,jj,kk] = ind2sub(size(ppgR.data),select_voxels);
ijk_func = [ii jj kk];
clear ii jj kk % housekeeping

% Get xyz coordinates of voxels in single subject space
xyz_anat = mrAnatXformCoords(acpcXform, ijk_func);

%%%%%%%%%%%% Render left and right hemispheres
hemis = {'r','l'};

for hh = 1:length(hemis)
    hemi_load = hemis{hh};
    
    gifti_name = fullfile(dDir,'derivatives','surfaces',['sub-' sub_label],['ses-' ses_label],['T1w_' hemi_load 'h_white_render.gii']);
    g = gifti(gifti_name);
    
    % get coordinates within only R hemipshere
    if ismember(ss,[2 3]) && isequal(hemi_load,'r')
        x_plot = -15; % midline off in this subject, should have acpc oriented...
        xyz_select = xyz_anat(:,1)>x_plot;
    elseif ~ismember(ss,[2 3]) && isequal(hemi_load,'r')
        x_plot = -10;
        xyz_select = xyz_anat(:,1)>x_plot;
    elseif ismember(ss,[2 3]) && isequal(hemi_load,'l')
        x_plot = 5; % midline off in this subject, should have acpc oriented...
        xyz_select = xyz_anat(:,1)<x_plot;
    elseif ~ismember(ss,[2 3]) && isequal(hemi_load,'l')
        x_plot = 10;
        xyz_select = xyz_anat(:,1)<x_plot;
    end
    xx_plot = xyz_anat(xyz_select,1);
    yy_plot = xyz_anat(xyz_select,2);
    zz_plot = xyz_anat(xyz_select,3);

    pc12_render_sel = pc12_render(xyz_select,:);
    
    maxPlot = 0.5;
    intensity_plot = pc12_render_sel./maxPlot;    
    intensity_plot(intensity_plot>1) = 1;
    intensity_plot(intensity_plot<-1) = -1;
    data_colors_rgb = bbData2Colors([intensity_plot(:,1) intensity_plot(:,2)]);
    
    figure
    brainHandle = bbRenderGifti(g); hold on
    for kk = 1:size(intensity_plot,1)
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',data_colors_rgb(kk,:))
    end
    % brainHandle.FaceAlpha = .5; % Make the brain transparent
    title(['R>' num2str(Rthreshold,3)])
    
    bbViewLight(90,0)
    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_view90_0']))
    
    bbViewLight(270,0)
    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_view270_0']))
    
    clear xyz_select g

end

%% Plot predicted responses for all subjects as function of time(s)
% this can be used as a colorscale to go with renderings

Rthreshold = .5;

% get a segmentation to create a brain, CSF and vessels mask:
% niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
niSegm2 = niftiRead([save_name_base '_combineSegm.nii.gz']);
brainMask = niSegm2.data>0;

% load PPG-RCOD
ppgR = niftiRead([save_name_base '_codPPG.nii.gz']);

% get model
load(fullfile(dDir,'derivatives','brainbeat','group','allsubs_pc12'),'pc1','pc2')
% get model timing
load(fullfile(dDir,'derivatives','brainbeat','group',['canonicalPC_leavout' sub_labels{ss}]),'t_svd')

% Use mask:
select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);
% PC1/PC2 weights for voxels to plot:
pc12_render = [pc1Weight.data(select_voxels) pc2Weight.data(select_voxels)];

pc1_plot = pc1;
pc2_plot = pc2;
% get heartrate to have interpretable timing in seconds again
ni = niftiRead(fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'func',...
    ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']));
physio      = physioCreate('nifti',ni);
ppg_cycle   = 1./physioGet(physio,'PPGrate');

% make time actual time in seconds
tt = t_svd*ppg_cycle;
% rotate to split into two groups
pc_complex = complex(pc12_render(:,1),pc12_render(:,2));
pc_angle = angle(pc_complex*(pi/10-1i)); % multiply by (1-i) to rotate 45 deg, pi/10-1i rotates a little further

% plot colorscale of all voxels rendered
maxPlot = 0.5;
intensity_plot = pc12_render./maxPlot;    
intensity_plot(intensity_plot>1) = 1;
intensity_plot(intensity_plot<-1) = -1;
data_colors_rgb = bbData2Colors([intensity_plot(:,1) intensity_plot(:,2)]);

figure('Position',[0 0 180 200])
for kk = 1:size(pc12_render,1)
    if pc_angle(kk)<0
        subplot(2,1,1),hold on
        x = pc12_render(kk,1);
        y = pc12_render(kk,2);
        plot(tt,x*pc1_plot + y*pc2_plot,'Color',data_colors_rgb(kk,:),'LineWidth',1)
    else
        subplot(2,1,2),hold on
        x = pc12_render(kk,1);
        y = pc12_render(kk,2);
        plot(tt,x*pc1_plot + y*pc2_plot,'Color',data_colors_rgb(kk,:),'LineWidth',1)
    end
end
set(gca,'FontName','Ariel')
subplot(2,1,1),xlim([-0.5 1.8])
subplot(2,1,2),xlim([-0.5 1.8])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',['subj' int2str(ss) '_run' int2str(rr) '_renderCM_predictedResp']))
print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','figures',['subj' int2str(ss) '_run' int2str(rr) '_renderCM_predictedResp']))



