clear all
close all

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
dDir = '/biac4/wandell/data/BrainBeat/data';

% chdir(dDir)

%% The T2* data are here.  

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

subj ='20141017_1242';    % Date _ Time out of NIMS
scan ='6_1_mux8fov4_r1_25s_4mmFA25';  % A particular data set
scan_code = '8202_6_1'; % code for saving this particular scan
fmri = fullfile(dDir,subj,scan, [scan_code '.nii.gz']);
ni = niftiRead(fmri);

%% Load the correlation with heartbeat (made with bbCorrelate2physio):

ppgRname = fullfile(dDir,subj,scan,[scan_code '_corrPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

%% Load the timeseries around PPG peak (made with bbResponse2physio):

ppgTSname = fullfile(dDir,subj,scan,[scan_code '_PPGtrigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scan_code '_PPGtrigResponseT.mat']),'t');

%% Get the anatomicals:

anat      ='9_1_T1w_1mm_sag';   % Anatomical data
anat      = fullfile(dDir,subj,anat,'8202_9_1.nii.gz');
niAnatomy = niftiRead(anat);

%%
%% coregister the functionals to the T1:
%%

%%%%% use the first nifti to allign, this one has the most structural info:
ni1=ni;
ni1.data=ni1.data(:,:,:,1);
ni1.dim=ni1.dim(1:3);
ni1.pixdim=ni1.pixdim(1:3);

niAnatomy.pixdim=niAnatomy.pixdim(1:3); % only uses the first three

% allign functionals to the T1:
acpcXform = dtiRawAlignToT1(ni1,niAnatomy,[], [], [], 1); % last 1 adds a figure
% this saves the reallignment matrix in the folder of the functionals, 

% if already done, load matrix:
% load(fullfile(dDir,subj,scan,[scan_code 'AcpcXform.mat']))


%% now overlay some plots

% then use dtiGetSlice to get the same slice from 2 sets
curPos = [1,10,-20]; 
sliceThisDim = 1; 
% imDims=[-90 -120 -60; 90 130 90];
imDims=[-90 -120 -120; 90 130 90];

% get a slice from the functionals
imgVol = ppgR.data;%ni.data(:,:,:,1);
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [4 4 4];
[imgSlice1,x1,y1,z1]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);

% get a slice from the timeSeries
imgVol = ppgTS.data;%ni.data(:,:,:,1);
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [4 4 4];
[imgSlice2,x2,y2,z2]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);

% get the same slice for the anatomy
imgVol = niAnatomy.data;
img2std = niAnatomy.qto_xyz;
sliceNum =curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [1 1 1];
[imgSlice,x,y,z]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);

% for x and y for plotting:
if sliceThisDim==1
    x1=z1; x=z;
elseif sliceThisDim==2
    y1=x1; y=x;
    x1=z1; x=z;
elseif sliceThisDim==3
    % x and y stay the same
end

% now plot stuff:
f=figure;
cm=colormap(jet);
close(f)

figure('Position',[0 0 1200 600])
subplot(1,2,1)
imagesc(x1(1,:),y1(:,1),imgSlice1)
hold on
axis image

subplot(1,2,2)
% show the background:
image(x(1,:),y(:,1),cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:)));
hold on
% colormap gray
axis image

% tranform imgSlice1 (r-values) to colormap values that I want to use
imgSlice1_color=cat(3,zeros(size(imgSlice1)),zeros(size(imgSlice1)),zeros(size(imgSlice1)));
for k_x = 1:length(x1(1,:))
    for k_y = 1:length(y1(:,1))
        imgSlice1_color(k_y,k_x,:) = cm(1+floor(abs(imgSlice1(k_y,k_x))*64),:); 
    end
end
% show the colormap
h=image(x1(1,:),y1(:,1),imgSlice1_color);
% and make them transparent
set(h,'AlphaData',.2*ones(size(imgSlice1)))

%% overlay with time series
f=figure;
cm=colormap(jet);
close(f)

figure('Position',[0 0 1200 600])

subplot(1,2,1)
% show the background:
imagesc(x1(1,:),y1(:,1),imgSlice1);
hold on
colormap gray
axis image

subplot(1,2,2)
% show the background:
image(x(1,:),y(:,1),cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:)));
hold on
% colormap gray
axis image

a=squeeze(imgSlice2(:,:,t>-.2 & t<1)); % 3rd dimension is time
r2_scale=sqrt(imgSlice1.^2);

for sub_p=1:2
    subplot(1,2,sub_p)
for k=1:size(a,1)
    for m=1:size(a,2)
        
        % get plotting location
        k_x=y1(k,1);
        m_y=x1(1,m);
        
        % subtract the mean:
        r_curve=squeeze(a(k,m,:)-mean(a(k,m,:)));
        
        % scale by strength of correlation with heartbeat:
%         r_curve=r_curve/max(abs(r_curve(:)));
%         r_curve=r_curve*r2_scale(k,m);
%         plot(k+[0:size(a,3)-1]/size(a,3)-.5,m+r_curve,'r')

        % here we add a weird scaling factor, to make the curve fit in the voxel size
        s_f=3; 
        
        % just scale by the maximum and color with correlation strength
        if max(r_curve)>.8
            r_curve=s_f*r_curve/max(r_curve);
        end
        c_use=cm(1+floor(r2_scale(k,m)*64),:);
        plot(m_y+s_f*[0:size(a,3)-1]/size(a,3)-2,k_x+r_curve,'Color',c_use)
        
        plot(m_y-2,k_x,'.','Color',[.5 .5 .5])
    end
end
end
