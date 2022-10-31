
%
% This script makes Figure 7 and Figure 8 from:
%
%%%%% Measuring brain beats: cardiac-aligned fast fMRI signals %%%%%
% Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell
%
%

clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%%

%% Plot Slope/Time/HWHM back on brain
%%
%%

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,1,1,1,1};

% get model for traces figures
load(fullfile(dDir,'derivatives','brainbeat','group','allsubs_pc12'),'pc1','pc2')

for ss = 4%1:6
    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);

    % Get base name of saved data
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);

    % get model timing
    load(fullfile(dDir,'derivatives','brainbeat','group',['canonicalPC_leavout' sub_labels{ss}]),'t_svd')
    
    % load ppgR, COD 
    % note the qto_xyz is not in T1 space, but volume matches
    ppgR = niftiRead([save_name_base '_codPPG.nii.gz']);

    % these are all in T1w space
    % load pc1 weights
    pc1Weight = niftiRead([save_name_base '_space-T1w_canoPc1Weights.nii']);
    % load pc2 weights
    pc2Weight = niftiRead([save_name_base '_space-T1w_canoPc2Weights.nii']);
    % load Slope
    svd_slope = niftiRead([save_name_base '_space-T1w_modelSlope.nii']);
    % load PeakTime
    svd_peakt = niftiRead([save_name_base '_space-T1w_modelOnset.nii']);
    % load FWHM
    svd_fwhm = niftiRead([save_name_base '_space-T1w_model_FWHM.nii']); 
    acpcXform = pc1Weight.qto_xyz;

    % Get functional voxels to plot with rendering
    Rthreshold = .5;

    % get a segmentation to create a brain, CSF and vessels mask:
    % niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
    niSegm2 = niftiRead([save_name_base '_combineSegm.nii.gz']);
    brainMask = niSegm2.data>0;

    % Use mask:
    select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);
    % Time, R and Slope for voxels to plot:
    ColorInt_render = [svd_peakt.data(select_voxels) ppgR.data(select_voxels) svd_slope.data(select_voxels) pc1Weight.data(select_voxels) pc2Weight.data(select_voxels)];

    % Get indiced of selected voxels
    [ii,jj,kk] = ind2sub(size(ppgR.data),select_voxels);
    ijk_func = [ii jj kk]; clear ii jj kk % housekeeping
    xyz_anat = mrAnatXformCoords(acpcXform, ijk_func); % xyz coordinates in subject space

    %%%%%%%%%%%% Render left hemisphere
    hemi_load = 'l';
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

    %%%% get predicted timeseries of all to be rendered voxels
    these_waveforms = ColorInt_render_sel(:,4:5)*[pc1 pc2]';

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

    % get heartrate to have interpretable timing in seconds again
    ni = niftiRead(fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']));
    physio      = physioCreate('nifti',ni);
    ppg_cycle   = 1./physioGet(physio,'PPGrate');
    % get actual time in seconds
    tt = t_svd*ppg_cycle;
    
    % get minmum positions for plotting and keeping scale constant
    xyz_min = min([xx_plot yy_plot zz_plot]);
    xyz_max = max([xx_plot yy_plot zz_plot]);

    for view_nrs = 1%1:2
        if view_nrs==1
            videoName = [save_name_base '_gifMesialView'];
            view_angle = [90 0];
        elseif view_nrs==2
            videoName = [save_name_base '_gifLateralView'];
            view_angle = [270 0];
        end
        fid = figure('Position',[0 0 600 335],'Color','w'); % correct size for twitter gif
     
    
        for this_t = 1:length(tt)%find(tt>(ppg_cycle/2),1) % just plot 1 cycle, as in Figure 6 and 7 ti
            % plot negative responses
            subplot(5,7,[1:3 8:10 15:17 22:24])
            brainHandle = bbRenderGifti(g); hold on
            plot3(xyz_min(1),xyz_min(2)-10,xyz_min(3)-10,'.','MarkerSize',100,'Color',[1 1 1])
            plot3(xyz_max(1),xyz_max(2)+10,xyz_max(3)+10,'.','MarkerSize',100,'Color',[1 1 1])
            plot3(xyz_min(1),xyz_max(2)+10,xyz_min(3)-10,'.','MarkerSize',100,'Color',[1 1 1])
            plot3(xyz_min(1),xyz_min(2)-10,xyz_max(3)+10,'.','MarkerSize',100,'Color',[1 1 1])
            % only plot size etc if signal exceeds threshold at time t
            for kk = 1:size(intensity_plot,1)
                if these_waveforms(kk,this_t)<-.2 % only plot size etc if signal exceeds threshold at time t
                    if ColorInt_render_sel(kk,3)<0 % negative slope: local minimum
                        c_use = squeeze(cm2D_blood(100,ceil(color_plot(kk)*99.5)+100,:));
                        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',max(round(intensity_plot(kk)*100)-50,1),'Color',c_use)
                    end
                end
            end
            bbViewLight(view_angle(1),view_angle(2))
    
            % plot positive responses
            subplot(5,7,[5:7 12:14 19:21 26:28])
            brainHandle = bbRenderGifti(g); hold on
            plot3(xyz_min(1),xyz_min(2)-10,xyz_min(3)-10,'.','MarkerSize',100,'Color',[1 1 1])
            plot3(xyz_max(1),xyz_max(2)+10,xyz_max(3)+10,'.','MarkerSize',100,'Color',[1 1 1])
            plot3(xyz_min(1),xyz_max(2)+10,xyz_min(3)-10,'.','MarkerSize',100,'Color',[1 1 1])
            plot3(xyz_min(1),xyz_min(2)-10,xyz_max(3)+10,'.','MarkerSize',100,'Color',[1 1 1])
          for kk = 1:size(intensity_plot,1)
                if these_waveforms(kk,this_t)>.2 % only plot size etc if signal exceeds threshold at time t
                    if ColorInt_render_sel(kk,3)>0 % positive slope: local maximum
                        c_use = squeeze(cm2D_csf(100,ceil(color_plot(kk)*99.5)+100,:));
                        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',max(round(intensity_plot(kk)*100)-50,1),'Color',c_use)
                    end
                end
            end
            bbViewLight(view_angle(1),view_angle(2))
        
            if view_nrs==1
                title({'Time with respect to PPG peak:',[int2str(tt(this_t)*1000) ' ms']},'Position',[-25 -120 110],'FontSize',10,'Interpreter','latex')
            elseif view_nrs==2
                title({'Time with respect to PPG peak:',[int2str(tt(this_t)*1000) ' ms']},'Position',[-25 120 110],'FontSize',10,'Interpreter','latex')
            end

            % Add waveforms/colors
            for kk = 1:size(intensity_plot,1)
                if ColorInt_render_sel(kk,3)<0 % negative slope: local minimum
                    if these_waveforms(kk,this_t)<-.2 % only plot size etc if signal exceeds threshold at time t
                        subplot(5,7,29:31),hold on
                        c_use = squeeze(cm2D_blood(100,ceil(color_plot(kk)*99.5)+100,:));
                        x = ColorInt_render_sel(kk,4); %PC1
                        y = ColorInt_render_sel(kk,5); %PC2
                        plot(tt,x*pc1 + y*pc2,'Color',c_use)
                    end
                elseif ColorInt_render_sel(kk,3)>0 % positive slope: local maximum
                    if these_waveforms(kk,this_t)>.2 % only plot size etc if signal exceeds threshold at time t
                        subplot(5,7,33:35),hold on
                        c_use = squeeze(cm2D_csf(100,ceil(color_plot(kk)*99.5)+100,:));
                        x = ColorInt_render_sel(kk,4); %PC1
                        y = ColorInt_render_sel(kk,5); %PC2
                        plot(tt,x*pc1 + y*pc2,'Color',c_use)
                    end
                end
            end
            subplot(5,7,29:31),hold on
            plot([tt(this_t) tt(this_t)],[-1 1],'k','LineWidth',2),set(gca,'YTickLabel','')
            plot(tt,-0.2+zeros(size(tt)),'k','LineWidth',2)
            ylim([-1 1]),xlim([min(tt) max(tt)])
            title({'Areas with local minima (blood vessels)','color = timing, size = reliablity','Waveforms of areas (2 cycles)'},'FontSize',10,'Interpreter','latex','Position',[0.4554 1.1 0])
            xlabel('Time (s)','Interpreter','latex')

            subplot(5,7,33:35),hold on
            plot([tt(this_t) tt(this_t)],[-1 1],'k','LineWidth',2),set(gca,'YTickLabel','')
            plot(tt,0.2+zeros(size(tt)),'k','LineWidth',2)
            ylim([-1 1]),xlim([min(tt) max(tt)])
            title({'Areas with local maxima (CSF spaces)','color = timing, size = reliablity','Waveforms of areas (2 cycles)'},'FontSize',10,'Interpreter','latex','Position',[0.4554 1.1 0])
            xlabel('Time (s)','Interpreter','latex')
    
            nr_frames = 1;
            % Write each frame to the file
            for m = 1:nr_frames % write X frames: decides speed
                exportgraphics(gcf,'./local/brainbeats_Hermesetal2022_mesial.gif','Append',true);
            end
            
            clf

        end

    end

end


