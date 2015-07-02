clear all
close all

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
dDir = '/biac4/wandell/data/BrainBeat/data';

% chdir(dDir)

%% The T2* data are here.  

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

s_nr = 2;
s_info = bb_subs(s_nr);
subj=s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

% Get the MRVenogram:
niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii.gz']));
% load Veno rotation matrix:
xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));


%%
%% now do an SVD on the odd responses:
%%

% load PPG responses
scan_nr = 2;
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};

% nifti:
ni=niftiRead(fullfile(dDir,subj,scan,[scanName]));

load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT']),'t')

% load average of all odd heartbeats:
ppgTS=niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_odd']));

% load average of all odd heartbeats:
ppgTSeven=niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_even']));

% load coregistration matrix:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform.mat']))

% Load the correlation with heartbeat (made with bbCorrelate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_corrPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

% scale the time-series matrix by the correlation
% divide by max %%%% maybe think about z-scoring instead
ppgTS.data = ppgTS.data ./ repmat(max(abs(ppgTS.data),[],4),[1,1,1,size(ppgTS.data,4)]);
ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
% multiply by correlation size (absolute)
ppgTS.data = ppgTS.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTS.data,4)]));
ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));

% reshape for SVD
a = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);

%%%%% TEST LOW PASS FILTER START
% srate = 20;
% band = 2;
% Rp   = 1; Rs = 20; % third order Butterworth
% high_p =  band(1)*2/srate;
% delta = 0.001*2/srate;
% high_s = min(1-delta,high_p+0.1);
% [n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
% [bf_b,bf_a] = butter(n_band,wn_band,'low');
% for k=1:size(a,1)
%     a(k,:) = filtfilt(bf_b,bf_a,double(a(k,:)));
% end
%%%%% TEST END

a = a(:,t>=-.5 & t<=2);
t_sel = t(t>=-.5 & t<=2);
a = a-repmat(mean(a,2),1,size(a,2)); % subtract the mean
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

%%%% put output in structures:
% put 2 components weights in a matrix
out = [];
for k=1:2
    out(k).weights = reshape(v(:,k),[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);
end

% %%%%% MODEL WITH 2 COMPONENTS:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';
svdResults.model = reshape(pred,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3) size(ppgTS.data,4)]);

%%
%% Get Voxel Curves
%%


coordName = 'SupSag1';
getCoord = [-6 -47 9]; % superior saggital1
coordName = 'SupSag2';
getCoord = [-2 -66 -25]; % superior saggital2
coordName = 'SupSag3';
getCoord = [1 -69 -37]; % superior saggital3
coordName = 'SupSag4';
getCoord = [4 -72 -57]; % superior saggital4
coordName = 'SupSag5';
getCoord = [7 -66 -70]; % superior saggital5
coordName = 'SupSag6';
getCoord = [3 -72 -70]; % superior saggital6
% coordName = 'TransR1';
% getCoord = [31 -48 -76]; % transverse R1
% coordName = 'TransL1';
% getCoord = [-30 -49 -89]; % transverse L1
% coordName = 'Straight1';
% getCoord = [-1 -42 -68]; % straight sinus1
% getCoord = [18 26 -78]; % cavernous R1
% getCoord = [-18 26 -78]; % cavernous L1

coordName = 'CerAqua1';
getCoord = [-3 2 -52]; % cerebral aquaduct1

% max for scaling:
scale2max = max(abs([out(1).weights(:); out(2).weights(:)]));
scale2max = .02;

%%% ONLY WORKS FOR DIMENSION 3, OTHERS NOT CORRECTLY SHOW VOXEL LOCATION!!!
% use dtiGetSlice to get the same slice from 2 sets
sliceThisDim = 3; %

if s_nr==2
    imDims = [-90 -120 -120; 90 130 90]; 
elseif s_nr==3
    imDims = [-90 -120 -100; 90 130 110]; 
end

curPos = getCoord; 

% functionals to ACPC space:
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [4 4 4];

% PC1 weight
imgVol = out(1).weights;
[imgSlicePC1,x1,y1,z1]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
% PC2 weight
imgVol = out(2).weights;
[imgSlicePC2]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
% model
imgVol = svdResults.model;
[imgSliceModel]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
% data
imgVol = ppgTS.data;
[imgSliceTS]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
imgVol = ppgTSeven.data;
[imgSliceTSeven]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
% R
imgVol = ppgR.data;
[imgSliceR]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);

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
imgVol = niVeno.data;%ni.data(:,:,:,1);
img2std = xf_veno.acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [1 1 1];
[imgSliceV]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);

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
    imgSlicePC1=imrotate(imgSlicePC1,90);
    imgSlicePC2=imrotate(imgSlicePC2,90);
    imgSliceTS=imrotate(imgSliceTS,90);
    imgSliceTSeven=imrotate(imgSliceTSeven,90);
    imgSliceModel=imrotate(imgSliceModel,90);
    imgSliceR=imrotate(imgSliceR,90);
    imgSlice=imrotate(imgSlice,90);
    imgSliceV=imrotate(imgSliceV,90);
elseif sliceThisDim==2
    y1=z1; y=z;
    % rotate the images
    imgSlicePC1=imrotate(imgSlicePC1,90);
    imgSlicePC2=imrotate(imgSlicePC2,90);
    imgSliceTS=imrotate(imgSliceTS,90);
    imgSliceTSeven=imrotate(imgSliceTSeven,90);
    imgSliceModel=imrotate(imgSliceModel,90);
    imgSliceR=imrotate(imgSliceR,90);
    imgSlice=imrotate(imgSlice,90);
    imgSliceV=imrotate(imgSliceV,90);
elseif sliceThisDim==3
    % x and y stay the same
end


% figure('Position',[0 0 300 400])
figure('Position',[0 0 900 800])

subplot(3,2,[1 3])
% show the background:
imagesc(x,y,imgSlice/max(imgSlice(:)));
colormap gray, set(gca,'CLim',[0 .3]), hold on ,axis image
% for k=1:size(imgSlice1,1)
%     for m=1:size(imgSlice1,2)
%         % get plotting location
%         k_x = y1(k);
%         m_y = x1(m);
%         val_plot = [imgSlice1(k,m) imgSlice2(k,m)]./scale2max;
%         % remove larger values to allow different max scaling
%         val_plot(abs(val_plot)>1)=sign(val_plot(abs(val_plot)>1))*1;
%         data_colors_rgb = bbData2Colors(val_plot);
%         plot(m_y,k_x,'.','Color',data_colors_rgb,'MarkerSize',12)
%     end
% end

plot(getCoord(1),getCoord(2),'ro')

subplot(3,2,[2 4])
% show the background:
imagesc(x,y,imgSliceV/max(imgSliceV(:)));
colormap gray, set(gca,'CLim',[0 .3]), hold on ,axis image

plot(getCoord(1),getCoord(2),'ro')

subplot(3,2,[5 6]),hold on
[~,x_ind]=min(abs(x1-getCoord(1)));
[~,y_ind]=min(abs(y1-getCoord(2)));
plot(t,squeeze(imgSliceModel(x_ind,y_ind,:)),'k')

plot(t,squeeze(imgSliceTS(x_ind,y_ind,:)),'r')
plot(t,squeeze(imgSliceTSeven(x_ind,y_ind,:)),'Color',[1 .7 .7])

title(['r=' num2str(imgSliceR(x_ind,y_ind))])

% ylim([-0.1 0.1])

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/test_pca_curves/' subj '_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_' coordName])

