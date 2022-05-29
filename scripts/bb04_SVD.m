
%
% This script runs the SVD for all individual subjects for the analysis in:
%
%%%%% Measuring brain beats: cardiac-aligned fast fMRI signals %%%%%
% Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell
%
% There are some options for making figures to check the outputs that are
% not used for manuscript figures
%

clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%% Get all relevant scans and PPG triggered response

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,[1 2],[1 2],[1 2],[1 2]};

for ss = 1:length(sub_labels) % subjects/ses/acq
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    
    for rr = 1%:length(run_nrs{ss}) % runs
        
        sub_label = sub_labels{ss};
        ses_label = ses_labels{ss};
        acq_label = acq_labels{ss};
        run_nr = run_nrs{ss}(rr);
        
        % Get the anatomicals:
        niAnatomy = niftiRead(fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'anat',...
            ['sub-' sub_label '_ses-' ses_label '_T1w.nii.gz']));
        
        % Get functional for physioGet
        ni = niftiRead(fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'func',...
            ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']));
        
        % Get PPG triggered curves
        save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
            ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);
        
        ppgTSodd = niftiRead([save_name_base '_PPGtrigResponse_odd.nii.gz']); % ppg triggered time series
        ppgTSeven = niftiRead([save_name_base '_PPGtrigResponse_even.nii.gz']); % ppg triggered time series
        ppgT = load([save_name_base '_PPGtrigResponseT.mat'],'t');
        t = ppgT.t;
        
        % Load coregistration matrix (for the functionals):
        load([save_name_base '_AcpcXform_new.mat']);
        acpcXform = acpcXform_new; clear acpcXform_new
        
        %%%% COD 
        ppgRname = [save_name_base '_codPPG.nii.gz'];
        ppgR = niftiRead(ppgRname); % correlation with PPG
        
        %%%% Do the SVD
        
        % Set maximum of ppgTS to 1 for each voxel
        ppgTSodd.data = ppgTSodd.data ./ repmat(max(abs(ppgTSodd.data),[],4),[1,1,1,size(ppgTSodd.data,4)]);
        ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
        % Multiply by COD size (absolute)
        ppgTSodd.data = ppgTSodd.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSodd.data,4)]));
        ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));
        
        % Reshape to voxel X time:
        train_set = reshape(ppgTSodd.data,[numel(ppgTSodd.data(:,:,:,1)) length(t)]);
        
        % Select times to include in SVD: we want -0.5 to 1.5 heartbeat cycle
        physio      = physioCreate('nifti',ni);
        ppg_cycle   = 1./physioGet(physio,'PPGrate');
        train_set   = train_set(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);
        t_sel       = t(t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);
        
        % Do the SVD:
        meanTS = mean(train_set,2);
        a = train_set-repmat(meanTS,1,size(train_set,2)); % subtract the mean
        a(isnan(a)) = 0; % replace NaN by zero
        [u,s,v] = svd(a','econ');
        s = diag(s);
        
        % Get cumulative explained variance:
        var_explained = cumsum(s.^2) / sum(s.^2);
        
        % Resample the eigenvectors to heartrate time and save for later use
        t_hr = linspace(min(t_sel),max(t_sel),128);
        y1 = interp1(t_sel,double(u(:,1)),t_hr);
        y2 = interp1(t_sel,double(u(:,2)),t_hr);
        y3 = interp1(t_sel,double(u(:,3)),t_hr);
        
        % test model fit on even heartbeats
        all_pred_acc = zeros(length(s),1);
        % test data
        test_set = reshape(ppgTSeven.data,[numel(ppgTSeven.data(:,:,:,1)) length(t)]);
        test_set = test_set(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);
        
        for mm = 1:length(s)
            % prediction:
            kk = 1:mm;
            % COD for 1:mm components
            pred_temp = [u(:,kk)*diag(s(kk))*v(:,kk)']';
            all_pred_acc(mm) = calccod(pred_temp(:),test_set(:),1,0,0)./100;
        end
        clear pred_temp
        
        save([save_name_base '_pc12'],'y1','y2','y3','t_hr','var_explained','all_pred_acc')
        
        %%%%% save outputs in nifti structures in T1 space
        
        % spatial weights pc1, pc2, pc3 to nifti structures:
        % put 2 components weights in a matrix of size x*y*z
        out = [];
        for kk = 1:3
            out(kk).weights = reshape(v(:,kk),[size(ppgTSodd.data,1) size(ppgTSodd.data,2) size(ppgTSodd.data,3)]);
            % save pc1 spatial weight:
            ni_save = ni;
            ni_save.data = out(kk).weights;
            % make sure T1 space:
            ni_save.qto_xyz = acpcXform;
            ni_save.qto_ijk = inv(acpcXform);
            ni_save.sto_xyz = acpcXform;
            ni_save.sto_ijk = inv(acpcXform);
            ni_save.fname = fullfile([save_name_base '_space-T1w_pc' int2str(kk) 'weights.nii.gz']);
            niftiWrite(ni_save,[ni_save.fname])
            clear ni_save
        end
        
        %%%%% MODEL WITH 2 COMPONENTS:
        pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';
        svdResults.model = reshape(pred,[size(ppgTSodd.data,1) size(ppgTSodd.data,2) size(ppgTSodd.data,3) size(pred,2)]);
        % test-retest error
        test_train_error = sqrt(sum((test_set - train_set).^2,2));
        % model error
        test_model_error = sqrt(sum((test_set - pred).^2,2));
        % relative RMS error:
        rel_rms_error = test_model_error./test_train_error;
        svdResults.error = reshape(rel_rms_error,[size(ppgTSodd.data,1) size(ppgTSodd.data,2) size(ppgTSodd.data,3)]);
        % save model
        ni_save = ni;
        ni_save.data = svdResults.model;
        ni_save.qto_xyz = acpcXform;
        ni_save.qto_ijk = inv(acpcXform);
        ni_save.sto_xyz = acpcXform;
        ni_save.sto_ijk = inv(acpcXform);
        ni_save.fname = fullfile([save_name_base '_space-T1w_modelpc12.nii.gz']);
        niftiWrite(ni_save,[ni_save.fname])
        % save error
        ni_save = ni;
        ni_save.data = svdResults.error;
        ni_save.qto_xyz = acpcXform;
        ni_save.qto_ijk = inv(acpcXform);
        ni_save.sto_xyz = acpcXform;
        ni_save.sto_ijk = inv(acpcXform);
        ni_save.fname = fullfile([save_name_base '_space-T1w_modelerrorpc12.nii.gz']);
        niftiWrite(ni_save,[ni_save.fname])
    end
end


%%
%% Now we make some figures to check things in one subject/run
%% these figures are not used in the manuscript


%% Plot 2 components:

L = length(u(:,1));
Fs = 1./mean(diff(t));
f = Fs * (0:(L/2))/L;

nr_pc_plot = 2;
pc_colors = {[0 .1 .5],[1 .5 0],'g','y'};
figure('Position',[0 0 260 270])

for kk = 1:nr_pc_plot

    % plot pc
    subplot(2,1,1),hold on
    plot(t_sel,u(:,kk),'Color',pc_colors{kk},'LineWidth',2)
    xlim([t_sel(1) t_sel(end)])
    xlabel('time (s)')
    title([acq_labels{ss}])

    % plot fft of pc
    p1 = abs(fft(u(:,kk))/L);
    p1 = p1(1:floor(L/2+1));
    p1(2:end-1) = 2*p1(2:end-1);
    
    subplot(2,2,3),hold on
    plot(f,p1,'Color',pc_colors{kk},'LineWidth',2)
    ylabel('|P(f)|')
    xlabel('frequency (Hz)')
       
end

subplot(2,1,1),hold on
% sum of 2 components squared
% plot(t_sel,sqrt(u(:,1).^2+u(:,2).^2),':','Color',[.5 .5 .5],'LineWidth',1)
set(gca,'XTick',[0 1])
subplot(2,2,3),hold on
legend({'pc1','pc2'})%,'pc3','pc4'})

%%% how many components do we want
subplot(2,2,4)
plot(all_pred_acc,'k')
title(['COD (R)'])
xlabel('nr of components')
ylim([0 1])
box off

%% plot pc1 weights versus pc2 weights

weight1 = s(1)*v(:,1);
weight2 = s(2)*v(:,2);
weight3 = s(3)*v(:,3);

weight1(ppgR.data<0.8) = [];
weight2(ppgR.data<0.8) = [];
weight3(ppgR.data<0.8) = [];

figure
plot3(weight1,weight2,weight3,'.')
axis equal
axis square
% color points by local density


%%
%% check: plot a number of components (spatial/temporal weights)
%%
nrc_plot = 7;

sl_plotx = 26;
sl_plotz = 20;

figure('Position',[0 0 1200 500])
for k=1:nrc_plot
    subplot(3,nrc_plot,k),hold on
%     plot(t_sel,ppgCurve/max(ppgCurve),'Color',[.5 .5 .5])
    plot(t_sel,u(:,k),'b')
    xlim([min(t_sel) max(t_sel)])
    title(['c' int2str(k) ' cumvar ' num2str(var_explained(k),2)])

    whole_brain_v = reshape(v(:,k),[size(ppgTSodd.data,1) size(ppgTSodd.data,2) size(ppgTSodd.data,3)]);

    subplot(3,nrc_plot,nrc_plot+k)
    imagesc(squeeze(whole_brain_v(:,:,sl_plotz)),[-.03 .03])
    axis image
    
    subplot(3,nrc_plot,2*nrc_plot+k)
    im_plot=imrotate(squeeze(whole_brain_v(sl_plotx,:,:)),90);
    imagesc(im_plot,[-.02 .02])
    clear im_plot
    axis image
end


%% Plot one component on anatomy 
niOverlay = ni; % just get the nifti structure of the functionals
w_plot = 1; % pc component
niOverlay.data = out(w_plot).weights;
sliceThisDim = 3;

if ss == 1
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [-10 50 -21];
    curPos = [-10 -20 -21];
    curPos = [1 -20 -21];
elseif ss == 2
    imDims = [-90 -120 -100; 90 130 110];
%     curPos = [1,4,38];
    curPos = [9,20,-12];
end

maxPlot = .01;
bbOverlayDotsAnat(niOverlay,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,maxPlot)
ylabel(['PC ' int2str(w_plot) ' max=' num2str(maxPlot,3)])


%% Plot shapes for PC1 and PC2 with 2D colorscale
maxPlot = .008;
t_select = (t>-0.2 & t<1);

%%%% Select voxels:
% Quick brain mask:
brain_vect = ni.data(:,:,:,4); brain_vect = brain_vect(:); 
% Correlation mask:
ppgR_vect = ppgR.data(:); 
select_voxels = brain_vect>10 & ppgR_vect>.6;

% Model prediction for selected voxels:
pred = [u(:,1:2)*diag(s(1:2))*v(select_voxels,1:2)']';

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
intensity_plot = v(select_voxels,1)./maxPlot;    
intensity_plot(intensity_plot>1) = 1;
intensity_plot(intensity_plot<-1) = -1;
color_plot = v(select_voxels,2)./maxPlot;
color_plot(color_plot>1) = 1;
color_plot(color_plot<-1) = -1;

figure('Position',[0 0 300 450]),hold on
for kk = 1:size(pred,1)
    if intensity_plot(kk)>0
        subplot(2,1,1),hold on
        c_use = squeeze(cm2D(ceil(intensity_plot(kk)*99+1),ceil(color_plot(kk)*99.5)+100,:));
        plot(t(t_select),pred(kk,t_select),'Color',c_use)
    elseif intensity_plot(kk)<0
        subplot(2,1,2),hold on
        c_use = squeeze(cm2D(ceil(-intensity_plot(kk)*99+1),ceil(-color_plot(kk)*99.5)+100,:));
        plot(t(t_select),pred(kk,t_select),'Color',c_use)       
    end
end
subplot(2,1,1), axis tight
title('PC1>0')
xlabel('Time (s)'),ylabel('Model prediction')
subplot(2,1,2), axis tight
title('PC1<0')
xlabel('Time (s)'),ylabel('Model prediction')
set(gcf,'PaperPositionMode','auto')



