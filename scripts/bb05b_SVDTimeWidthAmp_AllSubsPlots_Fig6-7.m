
%%
%% FIGURES & RENDERINGS in individual subjects
%%
%% Plot Slope/Time/HWHM back on brain
%%
%%

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

% get model for traces figures
load(fullfile(dDir,'derivatives','brainbeat','group','allsubs_pc12'),'pc1','pc2')
% get model timing
load(fullfile(dDir,'derivatives','brainbeat','group',['canonicalPC_leavout' sub_labels{ss}]),'t_svd')

for ss = 1:6;
    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);

    % Get base name of saved data
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

    % load ppgR, COD
    ppgR = niftiRead([save_name_base '_codPPG.nii.gz']);

    % load pc1 weights
    pc1Weight = niftiRead([save_name_base '_space-T1w_canoPc1Weights.nii.gz']);
    % load pc2 weights
    pc2Weight = niftiRead([save_name_base '_space-T1w_canoPc2Weights.nii.gz']);
    % load Slope
    svd_slope = niftiRead([save_name_base '_space-T1w_modelSlope.nii.gz']);
    % load PeakTime
    svd_peakt = niftiRead([save_name_base '_space-T1w_modelOnset.nii.gz']);
    % load FWHM
    svd_fwhm = niftiRead([save_name_base '_space-T1w_model_FWHM.nii.gz']); 

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

    % plot peakTime for slope<0
    acpcXform = svd_slope.qto_xyz;

    % show sagittal
    sliceThisDim = 1;

    % not implemented yet with current colorscale
    % % local minima
    % svd_peakt_plot =  svd_peakt; % time for color 
    % svd_intensity_plot = ppgR; % R for intensity
    % svd_intensity_plot.data(svd_slope.data>0) = 0; % remove positive
    % bbOverlayDotsAnat_Color2D(svd_peakt_plot,svd_intensity_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .1],12,0.00001)
    % set(gcf,'PaperPositionMode','auto')
    % print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_Neg_View' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))
    % 
    % % local maxima
    % svd_peakt_plot =  svd_peakt; % time for color 
    % svd_intensity_plot = ppgR; % R for intensity
    % svd_intensity_plot.data(svd_slope.data<0) = 0; % remove negative
    % bbOverlayDotsAnat_Color2D(svd_peakt_plot,svd_intensity_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,[50 .1],12,0.00001)
    % set(gcf,'PaperPositionMode','auto')
    % print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_Pos_View' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))

    % Get functional voxels to plot with rendering

    Rthreshold = .5;

    % get a segmentation to create a brain, CSF and vessels mask:
    % niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
    niSegm2 = niftiRead([save_name_base '_combineSegm.nii.gz']);
    brainMask = niSegm2.data>0;

    % Use mask:
    select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);
    % Time,R and Slope for voxels to plot:
    ColorInt_render = [svd_peakt.data(select_voxels) ppgR.data(select_voxels) svd_slope.data(select_voxels) pc1Weight.data(select_voxels) pc2Weight.data(select_voxels)];

    % Get indiced of selected voxels
    [ii,jj,kk] = ind2sub(size(ppgR.data),select_voxels);
    ijk_func = [ii jj kk]; clear ii jj kk % housekeeping
    xyz_anat = mrAnatXformCoords(acpcXform, ijk_func); % xyz coordinates in subject space

    %%%%%%%%%%%% Render left hemisphere
    hemi_load = 'l';
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);
    gifti_name = fullfile(dDir,'derivatives','surfaces',['sub-' sub_label],['ses-' ses_label],['T1w_' hemi_load 'h_white_render.gii']);
    g = gifti(gifti_name);

    % get coordinates within only left hemipshere
    if ss==2 || ss==3
        x_plot = 5; % midline off in this subject, should have acpc oriented...
    else
        x_plot = 10;
    end
    xyz_select = xyz_anat(:,1)<x_plot; % only left hemisphere
    xx_plot = xyz_anat(xyz_select,1);
    yy_plot = xyz_anat(xyz_select,2);
    zz_plot = xyz_anat(xyz_select,3);

    % Time, R, Slope, PC1, PC2, xyz_select gets left hemisphere
    ColorInt_render_sel = ColorInt_render(xyz_select,:); 

    % Set maximum for dot colors:
    maxPlotC = 50;
    maxPlotI = Rthreshold; % always plot at max color

    % Get colors for selected voxels
    color_plot = ColorInt_render_sel(:,1)./maxPlotC;
    color_plot(color_plot>1) = 1;
    color_plot(color_plot<-1) = -1;
    intensity_plot = ColorInt_render_sel(:,2);%./maxPlotI; 

    % Make 2D colormaps for blood and csf: one to vary color, the other varies intensity
    cm_csf = customcolormap([0 .45 .55 1], [.7 1 .6;.5 1 .9;.6 .1 .6;1 .9 .8], 200);
    cm_blood = customcolormap([0 .45 .55 1], [.4 1 1;.4 .4 1; 1 .4 .4;1 1 .4], 200);
    gray_vect = .2*ones(200,3);
    cm2D_csf = zeros(100,size(cm_csf,1),3);
    cm2D_blood = zeros(100,size(cm_blood,1),3);
    for kk = 1:100
        cm2D_csf(kk,:,:) = cm_csf*kk/100 + gray_vect*(100-kk)/100;
        cm2D_blood(kk,:,:) = cm_blood*kk/100 + gray_vect*(100-kk)/100;
    end

    % plot positive responses
    figure
    brainHandle = bbRenderGifti(g); hold on
    % brainHandle.FaceAlpha = .5; % Make the brain transparent
    for kk = 1:size(intensity_plot,1)
        if ColorInt_render_sel(kk,3)>0 % positive slope: local maximum
            c_use = squeeze(cm2D_csf(100,ceil(color_plot(kk)*99.5)+100,:));
            plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',round(intensity_plot(kk)*100)-40,'Color',c_use)
        end
    end
    title(['R>' num2str(Rthreshold,3)])
    bbViewLight(90,0)

    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','sub_render',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewPosMed_time']))
    bbViewLight(270,0)
    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','sub_render',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewPosLat_time']))


    % plot negative responses
    figure
    brainHandle = bbRenderGifti(g); hold on
    for kk = 1:size(intensity_plot,1)
        if ColorInt_render_sel(kk,3)<0 % negative slope: local minimum
            c_use = squeeze(cm2D_blood(100,ceil(color_plot(kk)*99.5)+100,:));
            plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',round(intensity_plot(kk)*100)-40,'Color',c_use)
        end
    end
    title(['R>' num2str(Rthreshold,3)])
    bbViewLight(90,0)
    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','sub_render',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewNegMed_time']))
    bbViewLight(270,0)
    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','sub_render',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewNegLat_time']))

    %%%% plot predicted timeseries of all rendered voxels

    % get heartrate to have interpretable timing in seconds again
    ni = niftiRead(fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']));
    physio      = physioCreate('nifti',ni);
    ppg_cycle   = 1./physioGet(physio,'PPGrate');
    % make time actual time in seconds
    tt = t_svd*ppg_cycle;

    figure('Position',[0 0 180 200])

    for kk = 1:size(intensity_plot,1)
        if ColorInt_render_sel(kk,3)>0 % positive slope: local maximum
            subplot(2,1,1),hold on
            c_use = squeeze(cm2D_csf(100,ceil(color_plot(kk)*99.5)+100,:));
            x = ColorInt_render_sel(kk,4); %PC1
            y = ColorInt_render_sel(kk,5); %PC2
            plot(tt,x*pc1 + y*pc2,'Color',c_use)
        elseif ColorInt_render_sel(kk,3)<0 % negative slope: local minimum
            subplot(2,1,2),hold on
            c_use = squeeze(cm2D_blood(100,ceil(color_plot(kk)*99.5)+100,:));
            x = ColorInt_render_sel(kk,4); %PC1
            y = ColorInt_render_sel(kk,5); %PC2
            plot(tt,x*pc1 + y*pc2,'Color',c_use)
        end
    end

    set(gca,'FontName','Ariel')
    subplot(2,1,1),xlim([-0.5 1.8])
    subplot(2,1,2),xlim([-0.5 1.8])

    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','sub_render',['subj' int2str(ss) '_run' int2str(rr) '_predictedResp']))
    print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group','sub_render',['subj' int2str(ss) '_run' int2str(rr) '_predictedResp']))
    
end


%%
%%
%%
%% Rest is old code:
%%
%% plot model of arteries, veins, ventricles

% load segmentation
niSegm = niftiRead([save_name_base '_r_DKTatlas_aseg.nii.gz']);
dkt_table = readtable('dkt_areas.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dkt_table_surface = readtable('dkt_areas_surface.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
roiNames = dkt_table.label;
roiCodes = dkt_table.label_nr;

Rthreshold = 0.5;
aca_dkt = dkt_table_surface.DKT_nr(dkt_table_surface.ind_arterial == 1001);
aca_codes = roiCodes(ismember(dkt_table.DKT_nr,aca_dkt)); % find matching volume index from dkt_table
aca_voxels = ismember(niSegm.data,aca_codes);

mca_dkt = dkt_table_surface.DKT_nr(dkt_table_surface.ind_arterial == 2001);
mca_codes = roiCodes(ismember(dkt_table.DKT_nr,mca_dkt)); % find matching volume index from dkt_table
mca_voxels = ismember(niSegm.data,mca_codes);

pca_dkt = dkt_table_surface.DKT_nr(dkt_table_surface.ind_arterial == 3001);
pca_codes = roiCodes(ismember(dkt_table.DKT_nr,pca_dkt)); % find matching volume index from dkt_table
pca_voxels = ismember(niSegm.data,pca_codes);

maxPlot = 1;
intensity_plot = [pc1Weight.data(:) pc2Weight.data(:)]./maxPlot;
intensity_plot(intensity_plot>1) = 1;
intensity_plot(intensity_plot<-1) = -1;

figure('Position',[0 0 150 100]),hold on
% subplot(3,1,1),hold on
plot_these_voxels = find(aca_voxels(:)==1 & ppgR.data(:)>=Rthreshold & pc1Weight.data(:)<0);
mean_plot = mean(intensity_plot(plot_these_voxels,1))*pc1_plot + mean(intensity_plot(plot_these_voxels,2))*pc2_plot;
data_colors_rgb = bbData2Colors([mean(intensity_plot(plot_these_voxels,1)) mean(intensity_plot(plot_these_voxels,2))]);
plot(tt,mean_plot,'Color',data_colors_rgb,'LineWidth',1)
plot_these_voxels = find(mca_voxels(:)==1 & ppgR.data(:)>=Rthreshold & pc1Weight.data(:)<0);
mean_plot = mean(intensity_plot(plot_these_voxels,1))*pc1_plot + mean(intensity_plot(plot_these_voxels,2))*pc2_plot;
data_colors_rgb = bbData2Colors([mean(intensity_plot(plot_these_voxels,1)) mean(intensity_plot(plot_these_voxels,2))]);
plot(tt,mean_plot,'Color',data_colors_rgb,'LineWidth',1)
plot_these_voxels = find(pca_voxels(:)==1 & ppgR.data(:)>=Rthreshold & pc1Weight.data(:)<0);
mean_plot = mean(intensity_plot(plot_these_voxels,1))*pc1_plot + mean(intensity_plot(plot_these_voxels,2))*pc2_plot;
data_colors_rgb = bbData2Colors([mean(intensity_plot(plot_these_voxels,1)) mean(intensity_plot(plot_these_voxels,2))]);
plot(tt,mean_plot,'Color',data_colors_rgb,'LineWidth',1)
    
% subplot(3,1,2),hold on
plot_these_voxels = find((niSegm.data(:)==4 | niSegm.data(:)==43 ) & ppgR.data(:)>=Rthreshold & pc1Weight.data(:)>0);
mean_plot = mean(intensity_plot(plot_these_voxels,1))*pc1_plot + mean(intensity_plot(plot_these_voxels,2))*pc2_plot;
data_colors_rgb = bbData2Colors([mean(intensity_plot(plot_these_voxels,1)) mean(intensity_plot(plot_these_voxels,2))]);
plot(tt,mean_plot,'Color',data_colors_rgb,'LineWidth',1)

% subplot(3,1,3),hold on
plot_these_voxels = find(niSegm2.data(:)==5 & ppgR.data(:)>=Rthreshold & pc1Weight.data(:)<0);
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
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_predictedRespRoisR' int2str(Rthreshold*100)]))
print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_predictedRespRoisR' int2str(Rthreshold*100)])) 

