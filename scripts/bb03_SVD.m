clear all
close all

%% Base data directory 
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% chdir(dDir)

%% The T2* data are here.  

% TODO: Also try first half and second half reliability

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

s_nr = 6;
s_info = bb_subs(s_nr);
subj = s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

% % Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));
% % load Veno rotation matrix:
% xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));


%%
%% SVD on the odd responses:
%%

data_in = 'PPG';

% load PPG responses
scan_nr = 1;
scan = s_info.scan{scan_nr};
scanName = s_info.scanName{scan_nr};

% nifti:
ni = niftiRead(fullfile(dDir,subj,scan,[scanName '.nii.gz']));

load(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponseT']),'t')

% load average of all odd heartbeats:
ppgTS = niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_odd.nii.gz']));

% load average of all even heartbeats:
ppgTSeven = niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_even.nii.gz']));

% Load coregistration matrix:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

%%%% Scale the time-series matrix by the reliability
% Get odd/even corr/corr (made with bbCod/Correlate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_cod' data_in '.nii.gz']);
% ppgRname = fullfile(dDir,subj,scan,[scanName '_corr' data_in '.nii.gz']);
ppgR = niftiRead(ppgRname); % COD between even and odd heartbeats
% Set maximum of ppgTS to 1 for each voxel
ppgTS.data = ppgTS.data ./ repmat(max(abs(ppgTS.data),[],4),[1,1,1,size(ppgTS.data,4)]);
ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
% Multiply by correlation size (absolute)
ppgTS.data = ppgTS.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTS.data,4)]));
ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));

% Reshape to voxel X time:
a = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);

% Select times to include in SVD:
physio      = physioCreate('nifti',ni);
ppg_cycle   = 1./physioGet(physio,'PPGrate');
a = a(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);
t_sel = t(t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);

% Do the SVD:
meanTS = mean(a,2);
a = a-repmat(meanTS,1,size(a,2)); % subtract the mean
[u,s,v] = svd(a','econ');
s = diag(s);

% Get to cumulative explained variance:
var_explained = cumsum(s.^2) / sum(s.^2);
clear a a_test

% Plot 2 components:
L = length(u(:,1));
Fs = 1./mean(diff(t));
f = Fs * (0:(L/2))/L;

nr_pc_plot = 2;
pc_colors = {'r','b','g','y'};
% figure('Position',[0 0 300 400])
figure('Position',[0 0 130 250])

for kk = 1:nr_pc_plot

    % plot pc
    subplot(2,1,1),hold on
    plot(t_sel,u(:,kk),'Color',pc_colors{kk},'LineWidth',2)
    xlim([t_sel(1) t_sel(end)])
    xlabel('time (s)')
    title(['FlipA = ' int2str(s_info.scanFA{scan_nr})])

    % plot fft of pc
    p1 = abs(fft(u(:,kk))/L);
    p1 = p1(1:floor(L/2+1));
    p1(2:end-1) = 2*p1(2:end-1);
    
    subplot(2,1,2),hold on
    plot(f,p1,'Color',pc_colors{kk},'LineWidth',2)
    ylabel('|P(f)|')
    xlabel('frequency (Hz)')
    
    
end
subplot(2,1,1),hold on
plot(t_sel,sqrt(u(:,1).^2+u(:,2).^2),':','Color',[.5 .5 .5],'LineWidth',1)
set(gca,'XTick',[0 1])
subplot(2,1,2),hold on
legend({'pc1','pc2'})%,'pc3','pc4'})
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/svd/pc1Amp_pc2Time/s' int2str(s_nr) '_scan' int2str(scan_nr) '_pc_fft'])
% print('-painters','-r300','-depsc',[dDir './figures/svd/pc1Amp_pc2Time/s' int2str(s_nr) '_scan' int2str(scan_nr) '_pc_fft'])

%% Resample the eigenvectors and test whether they span a subspace
%% Save for later use

figure
subplot(2,1,1),hold on
plot(u(:,1:2))

% y1 = resample(double(u(:,1)),128,length(t_sel));
% y2 = resample(double(u(:,2)),128,length(t_sel));
t_hr = linspace(min(t_sel),max(t_sel),128);
y1 = interp1(t_sel,double(u(:,1)),t_hr);
y2 = interp1(t_sel,double(u(:,2)),t_hr);

subplot(2,1,2),hold on
plot([y1' y2'])

save(['./local/s-' int2str(s_nr) '_scan-' int2str(scan_nr) 'pc12'],'y1','y2','t_hr')


%%
%% plot a number of components:
%%
nrc_plot=7;

sl_plotx = 26;
sl_plotz = 20;

figure('Position',[0 0 1200 500])
for k=1:nrc_plot
    subplot(3,nrc_plot,k),hold on
%     plot(t_sel,ppgCurve/max(ppgCurve),'Color',[.5 .5 .5])
    plot(t_sel,u(:,k),'b')
    xlim([min(t_sel) max(t_sel)])
    title(['c' int2str(k) ' cumvar ' num2str(var_explained(k),2)])

    whole_brain_v=reshape(v(:,k),[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);

    subplot(3,nrc_plot,nrc_plot+k)
    imagesc(squeeze(whole_brain_v(:,:,sl_plotz)),[-.03 .03])
    axis image
    
    subplot(3,nrc_plot,2*nrc_plot+k)
    im_plot=imrotate(squeeze(whole_brain_v(sl_plotx,:,:)),90);
    imagesc(im_plot,[-.02 .02])
    clear im_plot
    axis image
end

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/2014_12_presentation1/' subj '_' scan '_examples_SVDcurves'])

%%
%% now see how many components we want
%%
nr_c_test = 20;
all_pred_acc = zeros(nr_c_test,1);

% test data
a_test = reshape(ppgTSeven.data,[numel(ppgTSeven.data(:,:,:,1)) length(t)]);
a_test = a_test(:,t>=-.5 & t<=2);

for mm = 1:nr_c_test
    % prediction:
    kk = 1:mm;
    
%     % squared correlation:
%     vect = u(:,kk);
%     vect_w = v(:,kk);
%     pred = [vect * vect_w']';
%     all_pred_acc(mm) = corr(a_test(:),pred(:)).^2;

    % COD:
    pred = [u(:,kk)*diag(s(kk))*v(:,kk)']';
    all_pred_acc(mm) = calccod(pred(:),a_test(:),1,0,0)./100;
    
end

figure
subplot(2,1,1)
plot(all_pred_acc)
title(['COD (R)'])
xlabel('nr of components')
ylim([0 1])

[~,nr_c_keep] = max(all_pred_acc);
subplot(2,1,2)
plot(t_sel,u(:,[1:nr_c_keep]))
legend('c1','c2','c3','c4')
xlabel('time (s)')
title('components')

% print('-painters','-r300','-dpng',['./figures/test_pca/svd_comp_' subj '_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr)])

%% Plot the mean:
% Note that ppgTS is rescaled according to cod(R)

meanSig = ppgTS;
meanSig.data = reshape(meanTS,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);

sliceThisDim = 1;

if s_nr == 1
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-1,10,-20];
elseif s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-1,10,-20];
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
elseif s_nr == 4
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
end
bbOverlayFuncAnat(meanSig,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)

%% put output in structures:
% put 2 components weights in a matrix
out = [];
for k=1:3
    out(k).weights = reshape(v(:,k),[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);
end

% %%%%% MODEL WITH 2 COMPONENTS:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';
svdResults.model = reshape(pred,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3) size(pred,2)]);

%%%%% MODEL ERROR
% train and test-sets
train_set = reshape(ppgTS.data(:,:,:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle),[prod(ppgTS.dim(1:3)) size(pred,2)]);
test_set = reshape(ppgTSeven.data(:,:,:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle),[prod(ppgTSeven.dim(1:3)) size(pred,2)]);
% test-retest error
test_train_error = sqrt(sum((test_set - train_set).^2,2));
% model error
test_model_error = sqrt(sum((test_set - pred).^2,2));
% relative RMS error:
rel_rms_error = test_model_error./test_train_error;

svdResults.error = reshape(rel_rms_error,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);

% save first component weight:
% ni_save = ni;
% ni_save.data = out(1).weights;
% ni_save.fname = fullfile(dDir,subj,scan,[scanName '_pc1_weights']);
% niftiWrite(ni_save,[ni_save.fname])

%% Plot one component on anatomy 
niOverlay = ni;
w_plot = 1;
niOverlay.data = out(w_plot).weights;
sliceThisDim = 3;

if s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [-10 50 -21];
    curPos = [-10 -20 -21];
    curPos = [1 -20 -21];
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
%     curPos = [1,4,38];
    curPos = [9,20,-12];
end

maxPlot = .01;
bbOverlayDotsAnat(niOverlay,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,maxPlot)
% bbOverlayDotsVeno(niOverlay,niVeno,acpcXform,xf_veno.acpcXform,sliceThisDim,imDims,curPos,maxPlot)
ylabel(['PC ' int2str(w_plot) ' max=' num2str(maxPlot,3)])

%% Plot with PC 1 for intensity, PC 2 for color

maxPlot = .002;
% maxPlot = .001;

sliceThisDim = 3;
if s_nr == 1
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [0 26 17]; 
elseif s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [-12 50 -21];
%     curPos = [-10 -20 -21]; % for figure set
    curPos = [-11 35 -77]; % Carotid
%     curPos = [-2 26 -63]; % Basilar
%     curPos = [-1 26 -21]; % SliceThisDim 1 Anterior Cerebral Artery, used in example
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
%     curPos = [0,4,38];
%     curPos = [0,4,38]; % for figure set
    curPos = [1 26 -63]; % x = 1 SliceThisDim 1 for Anterior Cerebral Artery
    curPos = [1 20 -23]; % x = 1 SliceThisDim 1 for Anterior Cerebral Artery
elseif s_nr == 4
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [-2 4 30];%[x x 38]
elseif s_nr == 5
    imDims = [-90 -120 -100; 90 130 120];
    curPos = [1 18 53];
elseif s_nr == 6
    imDims = [-90 -120 -100; 90 130 120];
    curPos = [7 19 3];
elseif s_nr == 7
    imDims = [-90 -120 -100; 90 130 120];
    curPos = [19 40 -29];
end

% Make two figures, for PC1>0, one PC1<0. The size of PC1 indicates the
% color intensity, the color indicates the PC2 phase: timing
niColor = ni;
niIntensity = ni;
niColor.data = out(2).weights;
niIntensity.data = out(1).weights;
niIntensity.data(niIntensity.data<0) = 0; % only plot positive
bbOverlayDotsAnat_Color2D(niColor,niIntensity,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,maxPlot)
title(['PC1>0, slice ' int2str(curPos(sliceThisDim))])
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/svd/pc1Amp_pc2Time/ACA_subj' subj '_scan' int2str(scan_nr) '_PC1pos_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

niColor = ni;
niIntensity = ni;
niColor.data = -out(2).weights; % note that we have to flip color, to always give timing effects the same color
niIntensity.data = -out(1).weights;
niIntensity.data(niIntensity.data<0) = 0; 
bbOverlayDotsAnat_Color2D(niColor,niIntensity,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,maxPlot)
title(['PC1<0, slice ' int2str(curPos(sliceThisDim))])
% print('-painters','-r300','-dpng',[dDir './figures/svd/pc1Amp_pc2Time/ACA_subj' subj '_scan' int2str(scan_nr) '_PC1neg_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

% niColor = ni;
% niIntensity = ni;
% niColor.data = out(2).weights; % note that we have to flip color, to always give timing effects the same color
% niColor.data(out(1).weights<0) = -niColor.data(out(1).weights<0);% flip negative intensities to maintain phase color;
% niIntensity.data = abs(out(1).weights);
% bbOverlayDotsAnat_Color2D(niColor,niIntensity,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,maxPlot)
% title(['PC1>0 & PC1<0, slice ' int2str(curPos(sliceThisDim))])
% print('-painters','-r300','-dpng',[dDir './figures/svd/pc1Amp_pc2Time/ACA_subj' subj '_scan' int2str(scan_nr) '_PC1all_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

%% Plot PC1 versus PC2 weights

% Quick brain mask:
brain_vect = ni.data(:,:,:,4);
brain_vect = brain_vect(:);

% Correlation mask:
ppgR_vect = ppgR.data(:);

% Select volxels with decent reliability:
select_voxels = brain_vect>10 & ppgR_vect>.9;
% select_voxels = brain_vect>10 & ppgR_vect>.1;

figure,hold on
% plot(v(select_voxels,1),v(select_voxels,2),'ko')
set_1 = select_voxels & v(:,1)>=0; % CSF
set_2 = select_voxels & v(:,1)<0 & v(:,2)>0; % veins
set_3 = select_voxels & v(:,1)<0 & v(:,2)<0; % arteries
plot(v(set_1,1),v(set_1,2),'go')
plot(v(set_2,1),v(set_2,2),'ro')
plot(v(set_3,1),v(set_3,2),'bo')
axis equal

sliceThisDim = 1;
if s_nr == 1
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [0 26 17]; 
elseif s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [-12 50 -21];
%     curPos = [-10 -20 -21]; % for figure set
%     curPos = [-11 34 -71]; % Carotid
%     curPos = [-2 26 -63]; % Basilar
    curPos = [1 26 -21]; % SliceThisDim 1 Anterior Cerebral Artery, used in example
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
%     curPos = [0,4,38];
%     curPos = [0,4,38]; % for figure set
    curPos = [01 26 -63]; % x = 1 SliceThisDim 1 for Anterior Cerebral Artery
elseif s_nr == 4
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [-5 4 30];%[x x 38]
end

niColor = [1 0 0;0 0 1; 0 1 0];
R_th = 0.5;% .0001; % This does not make sense now?
niIntensity = ni;
niIntensity.data = out(1).weights;
niIntensity.data(ppgR.data<R_th) = 0;
niIntensity.data(out(1).weights>0 & ppgR.data>=R_th) = 3;
niIntensity.data(out(1).weights<0 & out(2).weights>0 & ppgR.data>=R_th) = 1;
niIntensity.data(out(1).weights<0 & out(2).weights<0 & ppgR.data>=R_th) = 2;

bbOverlayDotsAnat_PickColor(niIntensity,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,niColor);


%% Convert PC2 weight colorscale to timing:

maxPlot = .008;
% t_select = (t>-0.2 & t<1); % s = 2, slower heartrate
t_select = (t>-0.2 & t<.6); % s = 3, faster heartrate

%%%% Select voxels:
% Quick brain mask:
brain_vect = ni.data(:,:,:,4); brain_vect = brain_vect(:); 
% Correlation mask:
ppgR_vect = ppgR.data(:); 
select_voxels = brain_vect>10 & ppgR_vect>.5;

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

% Get timing for selected voxels
voxel_peakT = NaN(size(pred,1),1);
temp_t = t(t_select);
for kk = 1:size(pred,1) % voxel loop
    if intensity_plot(kk)>0
        [~,m_ind] = max(pred(kk,t_select));
        voxel_peakT(kk) = temp_t(m_ind);
    elseif intensity_plot(kk)<0
        [~,m_ind] = min(pred(kk,t_select));
        voxel_peakT(kk) = temp_t(m_ind);   
    end
end
clear temp_t

figure('Position',[0 0 300 450]),hold on
for kk = 1:size(pred,1)
    if intensity_plot(kk)>0
        subplot(2,1,1),hold on
        c_use = squeeze(cm2D(ceil(intensity_plot(kk)*99+1),ceil(color_plot(kk)*99.5)+100,:));
        plot(voxel_peakT(kk),intensity_plot(kk),'.','Color',c_use,'MarkerSize',20)
    elseif intensity_plot(kk)<0
        subplot(2,1,2),hold on
        c_use = squeeze(cm2D(ceil(-intensity_plot(kk)*99+1),ceil(-color_plot(kk)*99.5)+100,:));
        plot(voxel_peakT(kk),intensity_plot(kk),'.','Color',c_use,'MarkerSize',20)
    end
end
subplot(2,1,1),hold on
xlim([min(voxel_peakT)-.05 max(voxel_peakT)+.05 ])
subplot(2,1,2),hold on
xlim([min(voxel_peakT)-.05 max(voxel_peakT)+.05 ])

% print('-painters','-r300','-dpng',[dDir './figures/svd/pc1Amp_pc2Time/ACA_subj' subj '_scan' int2str(scan_nr) '_Color2Time'])
% print('-painters','-r300','-depsc',[dDir './figures/svd/pc1Amp_pc2Time/ACA_subj' subj '_scan' int2str(scan_nr) '_Color2Time'])

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
% Time series for selected voxels:
% pred = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);
% pred = pred(select_voxels,:);

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

% print('-painters','-r300','-dpng',[dDir './figures/svd/pc1Amp_pc2Time/subj' int2str(s_nr) '_scan' int2str(scan_nr) '_PC1and2pred_R0_6'])
% print('-painters','-r300','-depsc',[dDir './figures/svd/pc1Amp_pc2Time/subj' int2str(s_nr) '_scan' int2str(scan_nr) '_PC1and2pred_R0_6'])
% print('-painters','-r300','-dpng',[dDir './figures/svd/pc1Amp_pc2Time/subj' int2str(s_nr) '_scan' int2str(scan_nr) '_PC1and2rawTS_R0_6'])
% print('-painters','-r300','-depsc',[dDir './figures/svd/pc1Amp_pc2Time/subj' int2str(s_nr) '_scan' int2str(scan_nr) '_PC1and2rawTS_R0_6'])


%%
%% reliability
%%

% model:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';

train_set = reshape(ppgTS.data,[prod(ppgTS.dim(1:3)) ppgTS.dim(4)]);
test_set = reshape(ppgTSeven.data,[prod(ppgTSeven.dim(1:3)) ppgTSeven.dim(4)]);

% error
test_train_error = sqrt(sum((test_set - train_set).^2,2));

% model error
test_model_error = sqrt(sum((test_set - pred).^2,2));

% relative RMS error:
rel_rms_error = test_model_error./test_train_error;

figure,hist(rel_rms_error,50)

%%
%% plot PC shapes and colors
%%

s_vect = sqrt(v(:,1).^2 + v(:,2).^2);

% model prediction for all voxels:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';

% quick brain mask:
brain_vect = ni.data(:,:,:,4);
brain_vect = brain_vect(:);

% correlation mask:
ppgR_vect = ppgR.data(:);

% select_voxels = s_vect>0.01 & brain_vect>10 & abs(ppgR_vect)>.3;
select_voxels = brain_vect>10 & ppgR_vect>.9;
% select_voxels = brain_vect>10 & ppgR_vect>.1;

% get colors:
colors_plot = v(:,1:2);
colors_plot = colors_plot / max(abs(colors_plot(:)));
colors_plot = bbData2Colors(colors_plot);

figure
plot(v(select_voxels,1),v(select_voxels,2),'ko')
axis equal

v_plot = v(select_voxels,:);
voxel_nr = find(select_voxels>0);
colors_plot = colors_plot(select_voxels,:);

t_select = (t>-0.2 & t<1);

figure('Position',[0 0 700 700]),hold on
plot([-0.02 0.02],[0 0],'k')
plot([0 0],[-0.02 0.02],'k')
for k=1:length(v_plot) 
%     plot(v_plot(k,1),v_plot(k,2),'.','Color',colors_plot(k,:))

    plot(v_plot(k,1)+t(t_select)/2000,v_plot(k,2)+pred(voxel_nr(k),t_select)/2000,'Color',colors_plot(k,:))

end
axis equal

set(gcf,'PaperPositionMode','auto')
%     print('-painters','-r300','-dpng',fullfile(dDir,subj,scan,[scanName '_BBcurves_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))
% print('-painters','-r300','-dpng',['./figures/test_pca/svd_' subj '_pc12_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_plotcolorshape90'])

figure('Position',[0 0 300 450])
for kk = 1:length(v_plot) 
    if v_plot(kk,1)>0
        subplot(2,1,1),hold on
        plot(t(t_select),pred(voxel_nr(kk),t_select),'Color',colors_plot(kk,:))
    elseif v_plot(kk,1)<0
        subplot(2,1,2),hold on
        plot(t(t_select),pred(voxel_nr(kk),t_select),'Color',colors_plot(kk,:))
    end
end
subplot(2,1,1), title('first pc>0')
subplot(2,1,2), title('first pc<0')
axis tight
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/test_pca/svd_' subj '_pc12_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_plotcolorshapeflat'])

%%

figure('Position',[0 0 200 200]),hold on
% example data and figure
data_in = [0 0;0 1;1 0;1 1;0 -1;-1 0;-1 -1;-1 1;1 -1];
data_in = [data_in;.2 * data_in; .4 * data_in; .6 * data_in; .8 * data_in];
data_colors_rgb = bbData2Colors(data_in);

for k=1:length(data_in)
%     plot(data_in(k,1),data_in(k,2),'k.','MarkerSize',40)
%     plot(data_in(k,1),data_in(k,2),'.','Color',data_colors_rgb(k,:),'MarkerSize',35)
    pred_thisval = data_in(k,1)*u(:,1) + data_in(k,2)*u(:,2);
    plot(data_in(k,1)+t(t_select)/3,data_in(k,2)+pred_thisval(t_select)/3,'Color',data_colors_rgb(k,:),...
        'LineWidth',2)
end

xlabel('pc 1 weight')
ylabel('pc 2 weight')

axis equal

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/svd/2Dcolorspace/subj' subj '_scan' int2str(scan_nr) '_SVD2comp_colors'])
% print('-painters','-r300','-depsc',[dDir './figures/svd/2Dcolorspace/subj' subj '_scan' int2str(scan_nr) '_SVD2comp_colors'])

%% Vary PC1 and PC2 combinations:
t_select = (t>-0.2 & t<1);

figure('Position',[0 0 300 300])
% pretend data
% data_in = [(ones(11,1)) (-1:.2:1)'];
data_in = [(ones(11,1)) (-1:.2:1)'];

subplot(2,1,1),hold on
data_colors_rgb = bbData2Colors(data_in);

for k=1:length(data_in)
    pred_curve = data_in(k,1)*u(:,1) + data_in(k,2)*u(:,2);
    plot(t(t_select),pred_curve(t_select),...
        'Color',data_colors_rgb(k,:),...
        'LineWidth',2)
end

subplot(2,1,2),hold on
data_in(:,1) = -data_in(:,1);
data_colors_rgb = bbData2Colors(data_in);

for k=1:length(data_in)
    pred_curve = data_in(k,1)*u(:,1) + data_in(k,2)*u(:,2);
    plot(t(t_select),pred_curve(t_select),...
        'Color',data_colors_rgb(k,:),...
        'LineWidth',2)
end

