clear all
close all

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
% dDir = '/biac4/wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% chdir(dDir)

%% The T2* data are here.  

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

s_nr = 2;
s_info = bb_subs(s_nr);
subj=s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii']));

% Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii.gz']));
% % load Veno rotation matrix:
% xf_veno=load(fullfile(dDir,subj,veno,[veno_code 'AcpcXform.mat']));


%%
%% now do an SVD on the odd responses:
%%

data_in = 'PPG';

% load PPG responses
scan_nr = 2;
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};

% nifti:
ni=niftiRead(fullfile(dDir,subj,scan,[scanName '.nii.gz']));

load(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponseT']),'t')

% load average of all odd heartbeats:
ppgTS=niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_odd.nii.gz']));

% load average of all odd heartbeats:
ppgTSeven=niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_even.nii.gz']));

% load coregistration matrix:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% scale the time-series matrix by the correlation

% Load the correlation with heartbeat (made with bbCorrelate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_cod' data_in '.nii.gz']);
ppgR = niftiRead(ppgRname); % COD between even and odd heartbeats

% Set maximum of ppgTS to 1 for each voxel
ppgTS.data = ppgTS.data ./ repmat(max(abs(ppgTS.data),[],4),[1,1,1,size(ppgTS.data,4)]);
ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
% Multiply by correlation size (absolute)
ppgTS.data = ppgTS.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTS.data,4)]));
ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));

% reshape to voxel X time
a = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);

% select times to include in SVD
if isequal(data_in,'PPG')
    a = a(:,t>=-.5 & t<=2);
    t_sel = t(t>=-.5 & t<=2);
elseif isequal(data_in,'RESP')
    a = a(:,t>=-.5 & t<=4);
    t_sel = t(t>=-.5 & t<=4);
end
meanTS = mean(a,2);
a = a-repmat(meanTS,1,size(a,2)); % subtract the mean
[u,s,v]=svd(a','econ');
s=diag(s);

%%%% test for the sign of the 2nd component (also check for 1st???)
% in the 2nd component, the first peak should be negative
[~,pm_i]=findpeaks(double(u(:,2)),'minpeakdistance',10);
[~,nm_i]=findpeaks(-double(u(:,2)),'minpeakdistance',10);
first_peak = min([pm_i; nm_i]);
if u(first_peak,2)<0
    % do nothing
elseif u(first_peak,2)>0
    % reverse sign pc and weights
    u(:,2)=-u(:,2); v(:,2) = -v(:,2);
end

% get to cumulative explained variance:
var_explained=cumsum(s.^2) / sum(s.^2)

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
    % COD with yes to gain and yes to mean subtraction
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

%% check whether the second PC only shifts the first in time:

% Get the amplitude spectrum, if the amplitude spectrum is the same, and
% only the phase differs, the components are very similar, the second only
% shifts the first one back and forth in time:
L = length(u(:,1));
Fs = 1./mean(diff(t));
f = Fs * (0:(L/2))/L;

nr_pc_plot = 2;
pc_colors = {'r','b','g','y'};
figure('Position',[0 0 200 400])

for kk = 1:nr_pc_plot

    % plot pc
    subplot(2,1,1),hold on
    plot(t,u(:,kk),'Color',pc_colors{kk},'LineWidth',2)
    xlim([t(1) t(end)])
    xlabel('time (s)')

    % plot fft of pc
    p1 = abs(fft(u(:,kk))/L);
    p1 = p1(1:floor(L/2+1));
    p1(2:end-1) = 2*p1(2:end-1);
    
    subplot(2,1,2),hold on
    plot(f,p1,'Color',pc_colors{kk},'LineWidth',2)
    ylabel('|P(f)|')
    xlabel('frequency (Hz)')
end
legend({'pc1','pc2'})%,'pc3','pc4'})
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/svd/subj' subj '_scan' int2str(scan_nr) '_svdComponents_fft'])


%% Plot the mean:
% Note that ppgTS is rescaled according to cod(R)

meanSig = ppgTS;
meanSig.data = reshape(meanTS,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);

sliceThisDim = 3;

if s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [1,10,-20];
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
end
bbOverlayFuncAnat(meanSig,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)

%% put output in structures:
% put 2 components weights in a matrix
out = [];
for k=1:2
    out(k).weights = reshape(v(:,k),[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);
end

% %%%%% MODEL WITH 2 COMPONENTS:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';
svdResults.model = reshape(pred,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3) size(ppgTS.data,4)]);

%%%%% MODEL ERROR
% train and test-sets
train_set = reshape(ppgTS.data,[prod(ppgTS.dim(1:3)) ppgTS.dim(4)]);
test_set = reshape(ppgTSeven.data,[prod(ppgTSeven.dim(1:3)) ppgTSeven.dim(4)]);
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

%% plot a component on anatomy 
niOverlay = ni;
niOverlay.data = out(1).weights;
sliceThisDim = 3;

if s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-10 50 -21];
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
end

maxPlot = .01;
bbOverlayDotsAnat(niOverlay,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,maxPlot)

%% let's test: only plotting pc 2 for pc1<0/pc1>0
niOverlay = ni;
niOverlay.data = out(3).weights;
niOverlay.data(out(1).weights>0) = 0;
sliceThisDim = 1;
if s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-10 50 -21];
    curPos = [-1 50 -21];
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
end

maxPlot = .01;
bbOverlayDotsAnat(niOverlay,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,maxPlot)

%% plot component 1 and 2 combination on T1
%%
%% TODO: plot model shape - curve

% max for scaling:
scale2max = max(abs([out(1).weights(:); out(2).weights(:)]));
scale2max = .01;

% use dtiGetSlice to get the same slice from 2 sets
sliceThisDim = 3; 

if s_nr==2
    imDims = [-90 -120 -120; 90 130 90]; 
    curPos = [-10 50 -21];
elseif s_nr==3
    imDims = [-90 -120 -100; 90 130 110]; 
    curPos = [7,4,30]; 
end

% functionals to ACPC space
% settings:
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [4 4 4];

% functionals to ACPC
imgVol = out(1).weights;
[imgSlice1,x1,y1,z1]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
imgVol = out(2).weights;
[imgSlice2]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
if sliceThisDim == 1 || sliceThisDim == 3
    x1=x1(1,:)';
    y1=y1(:,1);
elseif sliceThisDim == 2
    x1=x1(:,1);
    y1=y1(1,:)';
end
z1=z1(1,:)';

% Anatomy to ACPC
imgVol = niAnatomy.data;
img2std = niAnatomy.qto_xyz;
sliceNum =curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [1 1 1];
[imgSlice,x,y,z]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
if sliceThisDim == 1 || sliceThisDim == 3
    x=x(1,:)';
    y=y(:,1);
elseif sliceThisDim == 2
    x=x(:,1);
    y=y(1,:)';
end
z=z(1,:)';

%%%% get the MRV instead of T1
%%%% get a slice from the MRV
% imgVol = niVeno.data;%ni.data(:,:,:,1);
% img2std = xf_veno.acpcXform;
% sliceNum = curPos(sliceThisDim);
% interpType = 'n';
% mmPerVox = [1 1 1];
% % mmPerVox = [4 4 4];
% [imgSlice,x,y,z]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
% if sliceThisDim == 1 || sliceThisDim == 3
%     x=x(1,:)';
%     y=y(:,1);
% elseif sliceThisDim == 2
%     x=x(:,1);
%     y=y(1,:)';
% end
% z=z(1,:)';
%%%% get the MRV instead of T1

% for x and y for plotting:
if sliceThisDim==1
    x1=z1; x=z;
    % flip x and y
    x1_t=x1;
    x1=y1; 
    y1=x1_t;
    x_t=x;
    x=y; 
    y=x_t;
    % and for the images
    imgSlice1=imrotate(imgSlice1,90);
    imgSlice2=imrotate(imgSlice2,90);
    imgSlice=imrotate(imgSlice,90);
elseif sliceThisDim==2
    y1=z1; y=z;
    % rotate the images
    imgSlice1=imrotate(imgSlice1,90);
    imgSlice2=imrotate(imgSlice2,90);
    imgSlice=imrotate(imgSlice,90);
elseif sliceThisDim==3
    % x and y stay the same
end


% figure('Position',[0 0 300 400])
figure('Position',[0 0 600 800])

% show the background:
imagesc(x,y,imgSlice/max(imgSlice(:)));
colormap gray
set(gca,'CLim',[0 .3])
hold on
axis image

for k=1:size(imgSlice1,1)
    for m=1:size(imgSlice1,2)
        % get plotting location
        k_x = y1(k);
        m_y = x1(m);

        val_plot = [imgSlice1(k,m) imgSlice2(k,m)]./scale2max;
        % remove larger values to allow different max scaling
        val_plot(abs(val_plot)>1)=sign(val_plot(abs(val_plot)>1))*1;
        data_colors_rgb = bbData2Colors(val_plot);
        
        plot(m_y,k_x,'.','Color',data_colors_rgb,'MarkerSize',12)

    end
end

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/svd/subj' subj '_scan' int2str(scan_nr) '_SVD2comp_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])


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
select_voxels = brain_vect>10 & ppgR_vect.^2>.9;

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

figure('Position',[0 0 300 300]),hold on
for k=1:length(v_plot) 
    plot(t(t_select),pred(voxel_nr(k),t_select),'Color',colors_plot(k,:))
end
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
% print('-painters','-r300','-dpng',[dDir './figures/svd/subj' subj '_scan' int2str(scan_nr) '_SVD2comp_colors'])
% print('-painters','-r300','-depsc',[dDir './figures/svd/subj' subj '_scan' int2str(scan_nr) '_SVD2comp_colors'])

%% Vary PC1 and PC2 combinations:

figure('Position',[0 0 300 300])
% example data and figure
% data_in = [0 0;0 1;1 0;1 1;0 -1;-1 0;-1 -1;-1 1;1 -1];

subplot(2,1,1),hold on
data_in = [(ones(11,1)) (-1:.2:1)'];
data_colors_rgb = bbData2Colors(data_in);

for k=1:length(data_in)
    pred_curve = data_in(k,1)*u(:,1) + data_in(k,2)*u(:,2);
    plot(t(t_select),pred_curve(t_select),...
        'Color',data_colors_rgb(k,:),...
        'LineWidth',2)
end

subplot(2,1,2),hold on
data_in = [(-ones(11,1)) (-1:.2:1)'];
data_colors_rgb = bbData2Colors(data_in);

for k=1:length(data_in)
    pred_curve = data_in(k,1)*u(:,1) + data_in(k,2)*u(:,2);
    plot(t(t_select),pred_curve(t_select),...
        'Color',data_colors_rgb(k,:),...
        'LineWidth',2)
end


%%
figure('Position',[0 0 700 700]),hold on
% example data and figure
data_in = [0 0;0 1;1 0;1 1;0 -1;-1 0;-1 -1;-1 1;1 -1];
data_in = [data_in;.2 * data_in; .4 * data_in; .6 * data_in; .8 * data_in];
data_colors_rgb = bbData2Colors(data_in);

for k=1:length(data_in)
%     plot(data_in(k,1),data_in(k,2),'k.','MarkerSize',40)
%     plot(data_in(k,1),data_in(k,2),'.','Color',data_colors_rgb(k,:),'MarkerSize',35)
    pred_thisval = data_in(k,1)*u(:,1) + data_in(k,2)*u(:,2);
    plot(t(t_select),pred_thisval(t_select),'Color',data_colors_rgb(k,:),...
        'LineWidth',2)
end

axis equal

set(gcf,'PaperPositionMode','auto')
%     print('-painters','-r300','-dpng',fullfile(dDir,subj,scan,[scanName '_BBcurves_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))
% print('-painters','-r300','-dpng',['./figures/test_pca/svd_' subj '_pc12_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_exampleshapes'])
