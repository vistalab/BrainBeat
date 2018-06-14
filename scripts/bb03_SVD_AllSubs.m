
% 1: We take the principle components, resampled to the heartrate for all
% subjects.
% 2: Do a PCA to see how they replicate across subjects.
% 3: Run through the data to test how well each voxel is predicted by a
% combination of these two components.
%
% DH & BW 2018, vistalab

%% Save canonical heartbeat responses

% subjects indices
s_nr = [2 3 4];
% scan nrs, here take the 25 or 20 degree flip angle
% scan_nr  = [1 1 2];
scan_nr  = [3 3 3];

% load first two principle components for these subjects and scans
all_pcs = zeros(length(s_nr),2,128);
for kk = 1:length(s_nr)
    load(['./local/s-' int2str(s_nr(kk)) '_scan-' int2str(scan_nr(kk)) 'pc12'],'y1','y2','t_hr')
    all_pcs(kk,1,:) = y1;
    all_pcs(kk,2,:) = y2;
end

t_svd = linspace(-.5,1.5,128);

% plot them:
figure
subplot(2,1,1),hold on
plot(t_svd,squeeze(all_pcs(:,1,:)),'k','LineWidth',2)
plot(t_svd,squeeze(all_pcs(:,2,:)),'Color',[1 .5 0],'LineWidth',2)
title('heartbeat components for 3 subjects')

% reshape into matrix size subject/scan*2 X time
all_pcs = reshape(all_pcs,size(all_pcs,1)*size(all_pcs,2),size(all_pcs,3))';

% do the pca on the first two pca's
[u,s,v] = svd(all_pcs);
% s = diag(s);

subplot(2,1,2),hold on
temp = u*s;
% plot(u(:,1:3))
plot(t_svd,temp(:,1),'k','LineWidth',2)
plot(t_svd,temp(:,2),'Color',[1 .5 0],'LineWidth',2)
xlabel('time (heartbeat cycles)')

pc1  = temp(:,1);
pc2  = temp(:,2);
pc3  = temp(:,3);
title('canonical heartbeat components')

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/svd/pc1Amp_pc2Time/canonicalPC_S234'])
print('-painters','-r300','-depsc',[dDir './figures/svd/pc1Amp_pc2Time/canonicalPC_S234'])

% save them:
save(['./local/allsubs_pc12'],'pc1','pc2','pc3')

%% Plot canonical PCs in figure:

load(['./local/allsubs_pc12'],'pc1','pc2','pc3')
pc1 = pc1(1:75);
pc2 = pc2(1:75);

figure('Position',[0 0 600 300])
subplot(2,2,1)
plot(pc1,'k','LineWidth',2)
axis tight
axis off
subplot(2,2,3)
plot(pc2,'k','LineWidth',2)
axis tight
axis off

% Make 2D colormap: one to vary color, the other varies intensity
cm = jet(250); cm = cm(26:225,:);
cm = cm(end:-1:1,:);
cm = cm+.4; cm(cm>1)=1;
gray_vect = .2*ones(200,3);
cm2D = zeros(100,size(cm,1),3);
for kk = 1:100
    cm2D(kk,:,:) = cm*kk/100 + gray_vect*(100-kk)/100;
end

subplot(1,2,2),hold on
for kk = -1:.4:1
    for ll = -1:.4:1
        if kk<0
            plot([kk:.3/74:kk+.3],ll+.3*(kk*pc1 + ll*pc2),...
                'Color',cm2D(round(-kk*100),round((-ll+1)*99+1),:),...
                'LineWidth',2)
        elseif kk>0
            plot([kk:.3/74:kk+.3],ll+.3*(kk*pc1 + ll*pc2),...
                'Color',cm2D(round(kk*100),round((ll+1)*99+1),:),...
                'LineWidth',2)
        end
    end
end
axis square
axis tight
axis off

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/svd/pc1Amp_pc2Time/model'])
print('-painters','-r300','-depsc',[dDir './figures/svd/pc1Amp_pc2Time/model'])


%% load canonical heartbeat responses and run through data
clear all
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% Load canonical heartbeat responses:
load(['./local/allsubs_pc12'],'pc1','pc2','pc3')

% Load PPG responses:
s_nr = 5;
scan_nr = 1;

s_info = bb_subs(s_nr);
subj=s_info.subj;
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};
data_in = 'PPG';

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

% Functionals:
ni=niftiRead(fullfile(dDir,subj,scan,[scanName '.nii.gz']));

% Load average of all odd heartbeats:
ppgTS=niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_odd.nii.gz']));

% Load average of all even heartbeats:
ppgTSeven=niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_even.nii.gz']));

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

% Load timing of heartbeat triggered responses:
load(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponseT']),'t')

% Load coregistration matrix:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% Times that were included in SVD:
physio      = physioCreate('nifti',ni);
ppg_cycle   = 1./physioGet(physio,'PPGrate');
t_sel = t(t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);

% Get even heartbeat responses that were included in SVD:
a = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);
a = a(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);
test_set = reshape(ppgTSeven.data,[prod(ppgTSeven.dim(1:3)) ppgTSeven.dim(4)]);
test_set = test_set(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);


% Resample canonical pc1 and pc2 to original timing:
t_hr = linspace(min(t_sel),max(t_sel),128);
y1 = interp1(t_hr,pc1,t_sel);
y2 = interp1(t_hr,pc2,t_sel);
% y3 = interp1(t_hr,pc3,t_sel);

disp('get beta weights')
beta_weights = zeros(size(a,1),2);
for kk = 1:size(a,1) % voxels
    if mod(kk,10000)==0, disp(['voxel ' int2str(kk) ' of ' int2str(size(a,1))]),end
    [B] = regress(a(kk,:)',[y1;y2]');
    beta_weights(kk,:) = B;
end
disp('done')

disp('get cod')
r_weights = zeros(size(a,1),1);
relRMS_weights = zeros(size(a,1),1);
% error
test_train_error = sqrt(sum((test_set - a).^2,2));
% model
model_v = beta_weights*[y1;y2];
% model error
test_model_error = sqrt(sum((test_set - model_v).^2,2));
% relative RMS error:
rel_rms_error = test_model_error./test_train_error;

for kk = 1:size(a,1) % voxels
    if mod(kk,10000)==0, disp(['voxel ' int2str(kk) ' of ' int2str(size(a,1))]),end
    % model_v = beta_weights(kk,:)*[y1;y2];
    r_weights(kk) = calccod(model_v(kk,:)',a(kk,:)',1,0,0)./100;
end
r_weights(isnan(r_weights)) = 0; % zero out NaN when model prediction is zeros
disp('done')

%% how good are canonical principle components

figure('Position',[0 0 200 150]),
[n,x] = hist(rel_rms_error,30);
bar(x,n,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0])
xlabel('relative root mean square error')
ylabel('number of voxels')
box off
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/renderCanonSvd/subj' int2str(s_nr) '_scan' int2str(scan_nr) '_relRMSE'])
print('-painters','-r300','-depsc',[dDir './figures/renderCanonSvd/subj' int2str(s_nr) '_scan' int2str(scan_nr) '_relRMSE'])


%%
voxel_nr = 28;
figure('Position',[0 0 200 150]),hold on
plot(t_sel,model_v(voxel_nr,:),'k','LineWidth',2)
plot(t_sel,a(voxel_nr,:),'b:','LineWidth',2)
plot(t_sel,test_set(voxel_nr,:),'r:','LineWidth',2)

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/renderCanonSvd/subj' int2str(s_nr) '_scan' int2str(scan_nr) '_pred_train_test'])
print('-painters','-r300','-depsc',[dDir './figures/renderCanonSvd/subj' int2str(s_nr) '_scan' int2str(scan_nr) '_pred_train_test'])


%% Plot R/Betas back on brain
niSeg = niftiRead(fullfile(dDir,subj,scan,[scanName '_combineSegm.nii.gz']));

% SPM segmentation
niSPM = niftiRead(fullfile(dDir,subj,scan,[scanName '_spmSeg.nii.gz']));
brain_vect = niSPM.data(:);

% Select voxels with decent reliability:
select_voxels = brain_vect>0 & r_weights>.7;
 
figure
subplot(2,1,1),hold on
plot(beta_weights(select_voxels,1),beta_weights(select_voxels,2),'ko')
% plot3(beta_weights(select_voxels,1),beta_weights(select_voxels,2),beta_weights(select_voxels,3),'ko')

% subplot(2,2,3),hold on
% set_1 = select_voxels & beta_weights(:,1)>=0; % CSF
% set_2 = select_voxels & beta_weights(:,1)<0 & beta_weights(:,2)>0; % veins
% set_3 = select_voxels & beta_weights(:,1)<0 & beta_weights(:,2)<0; % arteries
% plot(beta_weights(set_1,1),beta_weights(set_1,2),'go')
% plot(beta_weights(set_2,1),beta_weights(set_2,2),'ro')
% plot(beta_weights(set_3,1),beta_weights(set_3,2),'bo')
% axis equal
r_th = .5;
subplot(2,2,3),hold on
plot(beta_weights(niSeg.data(:)==5 & r_weights>r_th,1),beta_weights(niSeg.data(:)==5 & r_weights>r_th,2),'b.')
subplot(2,2,4),hold on
plot(beta_weights(niSeg.data(:)==3 & r_weights>r_th,1),beta_weights(niSeg.data(:)==3 & r_weights>r_th,2),'g.')
% plot(beta_weights(niSPM.data(:)==2 & r_weights>r_th,1),beta_weights(niSPM.data(:)==2 & r_weights>r_th,2),'ro')


sliceThisDim = 2;
if s_nr == 1
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [0 26 17]; 
elseif s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-1 50 -21]; % for figures 2
%     curPos = [-10 -20 -21]; % for figures 1
%     curPos = [-11 34 -71]; % Carotid
%     curPos = [-2 26 -63]; % Basilar
%     curPos = [1 26 -21]; % SliceThisDim 1 Anterior Cerebral Artery, used in example
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
%     curPos = [0,4,38];
%     curPos = [0,4,38]; % for figure set
    curPos = [1 26 -63]; % x = 1 SliceThisDim 1 for Anterior Cerebral Artery
elseif s_nr == 4
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0 4 35];%[x x 38] % x = 0 nicely captures posterior arteries/veins, x = -10, anterior middle artery
elseif s_nr == 5
    imDims = [-90 -120 -100; 90 130 120];
    curPos = [6,18,38];
end

ni_r = ni;
ni_r.data = ni_r.data(:,:,:,1);
ni_r.data(:) = r_weights;
ni_r.data(r_weights<.3) = 0;
ni_r.data(r_weights<0) = 0;
% bbOverlayDotsAnat(ni_r,niAnatomy,acpcXform,sliceThisDim,imDims,curPos);

ni_r = ppgR;
ni_r.data(ni_r.data<0.2) = 0;
% bbOverlayDotsAnat(ni_r,niAnatomy,acpcXform,sliceThisDim,imDims,curPos);

niColor = [1 0 0;0 0 1; 0 1 0];
R_th = 0.3;
niIntensity = ni;
niIntensity.data = zeros(size(ppgR.data));
niIntensity.data(ni_r.data<R_th) = 0;
niIntensity.data(beta_weights(:,1)>0 & r_weights>=R_th) = 3;
niIntensity.data(beta_weights(:,1)<0 & beta_weights(:,2)>0 & r_weights>=R_th) = 1;
niIntensity.data(beta_weights(:,1)<0 & beta_weights(:,2)<0 & r_weights>=R_th) = 2;


%%
figure
sliceThisDim = 1;
% curPos = [-20 -14 36];%[-1/-9 18 28] % zero
curPos = [6,18,38];
bbOverlayDotsAnat_PickColor(niIntensity,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,niColor);


%%
ni1 = ni;
ni1.data = ni1.data(:,:,:,1);
ni1.data(:) = beta_weights(:,1);

ni2 = ni;
ni2.data = ni2.data(:,:,:,1);
ni2.data(:) = beta_weights(:,2);
bbOverlayDotsAnat_FancyColorCircle(ni1,ni2,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.7);

%% plot beta weights for different tissue types

% Segmentation file:
segName = fullfile(dDir,subj,scan,[scanName '_combineSegm.nii.gz']);
niSeg = niftiRead(segName);

r_th = .5;
clear pc1w pc1w

for kk = 1:5
    pc1w(kk).beta = ni1.data(niSeg.data==kk & ppgR.data>r_th);
    pc2w(kk).beta = ni2.data(niSeg.data==kk & ppgR.data>r_th);
end

roiNames = {'GM','WM','Ventricles','CSF','Veno'};
figure('Position',[0 0 400 150])
for kk = 1:5
    subplot(1,5,kk),hold on
    histogram(pc1w(kk).beta,[-2.5:.4:2.5],'FaceColor',[.7 .7 .7])
    xlim([-2.5 2.5])
    title(roiNames{kk})
end
subplot(1,5,1)
ylabel(['number of voxels R>' num2str(r_th,3)])

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/svd/tissueWeights/PC1_subj' int2str(s_nr) '_scan' int2str(scan_nr) 'R' int2str(r_th*100)])
% print('-painters','-r300','-depsc',[dDir './figures/svd/tissueWeights/PC1_subj' int2str(s_nr) '_scan' int2str(scan_nr) 'R' int2str(r_th*100)])

%%
figure
for kk = 1:5
    subplot(1,5,kk),hold on
    hist(pc1w(kk).beta,[-2.5:.4:2.5])
    xlim([-2.5 2.5])
    title(roiNames{kk})
    
%     subplot(2,5,5+kk),hold on
%     hist(pc2w(kk).beta,[-2.5:.2:2.5])
% %     xlim([-2.5 2.5])
%     title(roiNames{kk})
end

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/reliable/MeanSig_subj' int2str(s_nr)])
% print('-painters','-r300','-depsc',[dDir './figures/reliable/MeanSig_subj' int2str(s_nr)])

%%
%% Get functional voxels to plot with rendering

Rthreshold = .5;

% SPM segmentation
niSPM = niftiRead(fullfile(dDir,subj,scan,[scanName '_spmSeg.nii.gz']));
brainMask = niSPM.data>0;

% Use mask:
select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);
% No mask:
% select_voxels = find(ppgR.data>=Rthreshold);

% Get indiced of selected voxels
[ii,jj,kk]=ind2sub(size(ppgR.data),select_voxels);
ijk_func=[ii jj kk];
clear ii jj kk % housekeeping

% Get xyz coordinates of voxels 
xyz_anat = mrAnatXformCoords(acpcXform, ijk_func);

% PC1/PC2 weights for voxels to plot:
pc12_render = [ni1.data(select_voxels) ni2.data(select_voxels)];

% Select hemisphere
xyz_select = xyz_anat(:,1)>-10;
xx_plot = xyz_anat(xyz_select,1);
yy_plot = xyz_anat(xyz_select,2);
zz_plot = xyz_anat(xyz_select,3);

pc12_render_sel = pc12_render(xyz_select,:);

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

%% Plot negative PC1 (veins/arteries) on render
hemi_load = 'r';
gifti_name = fullfile(dDir,subj,s_info.anat,['T1w_' hemi_load 'h_white_render.gii']);
% gifti_name = fullfile(dDir,subj,s_info.anat,['T1w_' hemi_load 'h_render.gii']);
g = gifti(gifti_name);

%%
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

% Just plot all points in red:
% plot3(xx_plot,yy_plot,zz_plot,'r.','MarkerSize',10)
bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/renderCanonSvd/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1neg_rh_lat'])

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/renderCanonSvd/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1neg_rh_med'])


% Plot positive PC1 (CSF)
figure
brainHandle=bbRenderGifti(g); hold on

for kk = 1:size(intensity_plot,1)
    if intensity_plot(kk)>0
        c_use = squeeze(cm2D(ceil(intensity_plot(kk)*99+1),ceil(color_plot(kk)*99.5)+100,:));
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',c_use)
    end
end

title(['R>' num2str(Rthreshold,3)])

% Just plot all points in red:
% plot3(xx_plot,yy_plot,zz_plot,'r.','MarkerSize',10)
bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/renderCanonSvd/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1pos_rh_lat'])

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/renderCanonSvd/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1pos_rh_med'])
