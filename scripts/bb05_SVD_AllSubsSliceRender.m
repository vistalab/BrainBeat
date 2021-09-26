
%%
%% FIGURES & RENDERINGS in individual subjects
%%
%% Plot R/Betas back on brain
%%
%%

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

% Get base name of saved data
save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

% load pc1
pc1Weight = niftiRead([save_name_base '_space-T1w_canoPc1Weights.nii.gz']);
% load pc2
pc2Weight = niftiRead([save_name_base '_space-T1w_canoPc2Weights.nii.gz']);
% load r weights (how much explained by canonical PC1-2)
r_weight = niftiRead([save_name_base '_space-T1w_canoPc12R.nii.gz']);

% Get anatomy
t1w_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'anat',...
            ['sub-' sub_label '_ses-' ses_label '_T1w.nii.gz']);
niAnatomy = niftiRead(fullfile(dDir,t1w_BIDSname));

if ss == 1
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-4 26 17]; 
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
    curPos = [-2,18,38]; % [-2,18,38]
elseif ss == 6
    imDims = [-90 -120 -100; 90 130 120];
    curPos = [-4,18,38];
end

% plot beta1 (pc1) and beta2 (pc2) using fancy color circle
acpcXform = pc1Weight.qto_xyz;

% plot entire circle
sliceThisDim = 1;
bbOverlayDotsAnat_FancyColorCircle(pc1Weight,pc2Weight,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.7);
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_exampleAxial38']))

% % now only select 7.30-1.30 on a clock
% % do this to test:
% % x = [1:-0.1:0 0:-0.1:-1 -1:0.1:0 0:0.1:1];
% % y = [0:0.1:1 1:-0.1:0 0:-0.1:-1 -1:0.1:0];
% % b = complex(x,y)
% % figure,plot(angle(b*(pi/4-pi/4*1i)))
% 
% pc_complex = complex(pc1Weight.data,pc2Weight.data);
% pc_angle = angle(pc_complex*(1-i)); % multiply by (1-i) to rotatio 45 deg
% 
% % plot veins and arteries
% pc1_plot = pc1Weight;
% pc2_plot = pc2Weight;
% pc1_plot.data(pc_angle<0) = 0;
% pc2_plot.data(pc_angle<0) = 0;
% bbOverlayDotsAnat_FancyColorCircle(pc1_plot,pc2_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.7);
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_example_pc12veins']))
% print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_example_pc12veins']))
% 
% % plot csf
% pc1_plot = pc1Weight;
% pc2_plot = pc2Weight;
% pc1_plot.data(pc_angle>0) = 0;
% pc2_plot.data(pc_angle>0) = 0;
% bbOverlayDotsAnat_FancyColorCircle(pc1_plot,pc2_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.7);
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_example_pc12csf']))
% print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_example_pc12csf']))

%%
%% Get functional voxels to plot with rendering

Rthreshold = .5;

% get a segmentation to create a brain, CSF and vessels mask:
% niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
niSegm2 = niftiRead([save_name_base '_combineSegm.nii.gz']);
brainMask = niSegm2.data>0;

% load PPG-RCOD
ppgR = niftiRead([save_name_base '_codPPG.nii.gz']);

% Use mask:
select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);
% PC1/PC2 weights for voxels to plot:
pc12_render = [pc1Weight.data(select_voxels) pc2Weight.data(select_voxels)];

% Get indiced of selected voxels
[ii,jj,kk] = ind2sub(size(ppgR.data),select_voxels);
ijk_func = [ii jj kk];
clear ii jj kk % housekeeping

% Get xyz coordinates of voxels in single subject space
xyz_anat = mrAnatXformCoords(acpcXform, ijk_func);


%%%%%%%%%%%% Render right hemisphere
hemi_load = 'r';
save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

gifti_name = fullfile(dDir,'derivatives','surfaces',['sub-' sub_label],['ses-' ses_label],['T1w_' hemi_load 'h_white_render.gii']);
g = gifti(gifti_name);

% get coordinates within only this hemipshere
% xyz_select = xyz_anat(:,1)<10;
if ss==2 || ss==3
    x_plot = -15; % midline off in this subject, should have acpc oriented...
else
    x_plot = -10;
end
xyz_select = xyz_anat(:,1)>x_plot;
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
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewLat']))

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewMed']))

%%%%%%%%%%%%% Render left hemisphere
hemi_load = 'l';
save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

gifti_name = fullfile(dDir,'derivatives','surfaces',['sub-' sub_label],['ses-' ses_label],['T1w_' hemi_load 'h_white_render.gii']);
g = gifti(gifti_name);

% get coordinates within only this hemipshere
if ss==2 || ss==3
    x_plot = 5; % midline off in this subject, should have acpc oriented...
else
    x_plot = 10;
end
xyz_select = xyz_anat(:,1)<x_plot;
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
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewMed']))

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewLat']))

%% Plot predicted responses for all subjects as function of time(s)
% this can be used as a colorscale to go with renderings

Rthreshold = .5;

% get a segmentation to create a brain, CSF and vessels mask:
% niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
niSegm2 = niftiRead([save_name_base '_combineSegm.nii.gz']);
brainMask = niSegm2.data>0;

% load PPG-RCOD
ppgR = niftiRead([save_name_base '_codPPG.nii.gz']);

% Use mask:
select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);
% PC1/PC2 weights for voxels to plot:
pc12_render = [pc1Weight.data(select_voxels) pc2Weight.data(select_voxels)];

maxPlot = 0.5;
intensity_plot = pc12_render./maxPlot;    
intensity_plot(intensity_plot>1) = 1;
intensity_plot(intensity_plot<-1) = -1;
data_colors_rgb = bbData2Colors([intensity_plot(:,1) intensity_plot(:,2)]);

% get model
load(fullfile(dDir,'derivatives','brainbeat','group','allsubs_pc12'),'pc1','pc2')
% get model timing
load(fullfile(dDir,'derivatives','brainbeat','group',['canonicalPC_leavout' sub_labels{ss}]),'t_svd')

pc1_plot = pc1;
pc2_plot = pc2;
% get heartrate to have interpretable timing in seconds again
ni = niftiRead(fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'func',...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']));
physio      = physioCreate('nifti',ni);
ppg_cycle   = 1./physioGet(physio,'PPGrate');

% make time actual time in seconds
tt = t_svd*ppg_cycle;
% rotate to split into two groups
pc_complex = complex(pc12_render(:,1),pc12_render(:,2));
pc_angle = angle(pc_complex*(pi/10-1i)); % multiply by (1-i) to rotate 45 deg, pi/10-1i rotates a little further

%% plot model of arteries, veins, ventricles

Rthreshold = .5;

% load segmentation
niSegm = niftiRead([save_name_base '_r_DKTatlas_aseg.nii.gz']);
dkt_table = readtable('dkt_areas.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dkt_table_surface = readtable('dkt_areas_surface.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
roiNames = dkt_table.label;
roiCodes = dkt_table.label_nr;

aca_dkt = dkt_table_surface.DKT_nr(dkt_table_surface.ind_arterial == 1001);
aca_codes = roiCodes(ismember(dkt_table.DKT_nr,aca_dkt)); % find matching volume index from dkt_table
aca_voxels = ismember(niSegm.data,aca_codes);

mca_dkt = dkt_table_surface.DKT_nr(dkt_table_surface.ind_arterial == 2001);
mca_codes = roiCodes(ismember(dkt_table.DKT_nr,mca_dkt)); % find matching volume index from dkt_table
mca_voxels = ismember(niSegm.data,mca_codes);

pca_dkt = dkt_table_surface.DKT_nr(dkt_table_surface.ind_arterial == 3001);
pca_codes = roiCodes(ismember(dkt_table.DKT_nr,pca_dkt)); % find matching volume index from dkt_table
pca_voxels = ismember(niSegm.data,pca_codes);

% pc12_plot = [pc1Weight.data(aca_voxels) pc2Weight.data(aca_voxels)];
maxPlot = 0.5;
intensity_plot = [pc1Weight.data(:) pc2Weight.data(:)]./maxPlot;
intensity_plot(intensity_plot>1) = 1;
intensity_plot(intensity_plot<-1) = -1;

figure('Position',[0 0 180 200]),hold on
% subplot(3,1,1),hold on
plot_these_voxels = find(aca_voxels(:)==1 & ppgR.data(:)>=Rthreshold);
mean_plot = mean(intensity_plot(plot_these_voxels,1))*pc1_plot + mean(intensity_plot(plot_these_voxels,2))*pc2_plot;
data_colors_rgb = bbData2Colors([mean(intensity_plot(plot_these_voxels,1)) mean(intensity_plot(plot_these_voxels,2))]);
plot(tt,mean_plot,'Color',data_colors_rgb,'LineWidth',1)
plot_these_voxels = find(mca_voxels(:)==1 & ppgR.data(:)>=Rthreshold);
mean_plot = mean(intensity_plot(plot_these_voxels,1))*pc1_plot + mean(intensity_plot(plot_these_voxels,2))*pc2_plot;
data_colors_rgb = bbData2Colors([mean(intensity_plot(plot_these_voxels,1)) mean(intensity_plot(plot_these_voxels,2))]);
plot(tt,mean_plot,'Color',data_colors_rgb,'LineWidth',1)
plot_these_voxels = find(pca_voxels(:)==1 & ppgR.data(:)>=Rthreshold);
mean_plot = mean(intensity_plot(plot_these_voxels,1))*pc1_plot + mean(intensity_plot(plot_these_voxels,2))*pc2_plot;
data_colors_rgb = bbData2Colors([mean(intensity_plot(plot_these_voxels,1)) mean(intensity_plot(plot_these_voxels,2))]);
plot(tt,mean_plot,'Color',data_colors_rgb,'LineWidth',1)
    
% subplot(3,1,2),hold on
plot_these_voxels = find((niSegm.data(:)==4 | niSegm.data(:)==43 ) & ppgR.data(:)>=Rthreshold);
mean_plot = mean(intensity_plot(plot_these_voxels,1))*pc1_plot + mean(intensity_plot(plot_these_voxels,2))*pc2_plot;
data_colors_rgb = bbData2Colors([mean(intensity_plot(plot_these_voxels,1)) mean(intensity_plot(plot_these_voxels,2))]);
plot(tt,mean_plot,'Color',data_colors_rgb,'LineWidth',1)

% subplot(3,1,3),hold on
plot_these_voxels = find(niSegm2.data(:)==5 & ppgR.data(:)>=Rthreshold);
mean_plot = mean(intensity_plot(plot_these_voxels,1))*pc1_plot + mean(intensity_plot(plot_these_voxels,2))*pc2_plot;
data_colors_rgb = bbData2Colors([mean(intensity_plot(plot_these_voxels,1)) mean(intensity_plot(plot_these_voxels,2))]);
plot(tt,mean_plot,':','Color',data_colors_rgb,'LineWidth',1)
% all_plot = bsxfun(@times,intensity_plot(plot_these_voxels,1)',pc1_plot) + bsxfun(@times,intensity_plot(plot_these_voxels,2)',pc2_plot);
% all_colors = bbData2Colors(intensity_plot(plot_these_voxels,:));
% for kk = 1:size(all_plot,2)
%     plot(tt,all_plot(:,kk),'Color',all_colors(kk,:),'LineWidth',1)
% end

set(gca,'FontName','Ariel')
% subplot(3,1,1),xlim([-0.5 1.8])
% subplot(3,1,2),xlim([-0.5 1.8])
% subplot(3,1,3),xlim([-0.5 1.8])

xlim([-0.5 1.8])
yline(0)
for kk = -0.5:.1:1.8
    plot([kk kk],[-0.02 0.02],'Color',[0 0 0])
end
plot([0 0],[-0.05 0.05],'Color',[0 0 0])
plot([.5 .5],[-0.05 0.05],'Color',[0 0 0])
plot([1 1],[-0.05 0.05],'Color',[0 0 0])
plot([1.5 1.5],[-0.05 0.05],'Color',[0 0 0])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_predictedRespRois']))
print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_predictedRespRois'])) 


%% %% plot colorscale of all voxels rendered
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
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_predictedResp']))
print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_predictedResp']))





%%
%%
%%
%%
%% only plot upper left/lower right quadrant
%%
pc_complex = complex(intensity_plot(:,1),intensity_plot(:,2));
pc_angle = angle(pc_complex*(1-1i)); % multiply by (1-i) to rotatio 45 deg
figure
brainHandle = bbRenderGifti(g); hold on
for kk = 1:size(intensity_plot,1)
    if pc_angle(kk)<0
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',data_colors_rgb(kk,:))
    end
end
title(['R>' num2str(Rthreshold,3)])


%%
%% plot with 2D red/blue colormap
%%

% Set maximum for dot colors:
maxPlot = .5;

% Make 2D colormap: one to vary color, the other varies intensity
cm = jet(250); cm = cm(26:225,:);
cm = cm(end:-1:1,:);
cm = cm+.4; cm(cm>1)=1;
gray_vect = .2*ones(200,3);
cm2D = zeros(100,size(cm,1),3);
for kk = 1:100
    cm2D(kk,:,:) = cm*kk/100 + gray_vect*(100-kk)/100;
end

% Get colors for selected voxels
intensity_plot = pc12_render_sel(:,1)./maxPlot;    
intensity_plot(intensity_plot>1) = 1;
intensity_plot(intensity_plot<-1) = -1;
color_plot = pc12_render_sel(:,2)./maxPlot;
color_plot(color_plot>1) = 1;
color_plot(color_plot<-1) = -1;

figure
brainHandle = bbRenderGifti(g); hold on
% brainHandle.FaceAlpha = .5; % Make the brain transparent

for kk = 1:size(intensity_plot,1)
    if intensity_plot(kk)<0
        c_use = squeeze(cm2D(ceil(-intensity_plot(kk)*99+1),ceil(-color_plot(kk)*99.5)+100,:));
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',c_use)
    end
end

title(['R>' num2str(Rthreshold,3)])

bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','renderCanonSvd',['sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1neg_rh_lat']))

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','renderCanonSvd',['sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1neg_rh_med']))

%% Plot positive PC1 (CSF)
figure
brainHandle = bbRenderGifti(g); hold on

for kk = 1:size(intensity_plot,1)
    if intensity_plot(kk)>0
        c_use = squeeze(cm2D(ceil(intensity_plot(kk)*99+1),ceil(color_plot(kk)*99.5)+100,:));
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',c_use)
    end
end

title(['R>' num2str(Rthreshold,3)])

bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','renderCanonSvd',['sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1pos_rh_lat']))

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','renderCanonSvd',['sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1pos_rh_med']))
