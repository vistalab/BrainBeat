function [] = bbOverlayFuncAnat(ni,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,varargin)
% Function to plot a functional and overlay with the anatomy

% Inputs: 
%   ni: functional scan
%   niAnatomy: anatomical scan
%   acpcXform: transformation matrix from nifti ijk to anatomy mm
%
%   varargin{1} = nr of figures: 2 plots the overlay additionally in a separate figure
%
%%%% now overlay r-map with anatomy

% % then use dtiGetSlice to get the same slice from 2 sets
% sliceThisDim = 2;
% if s==2
%     imDims = [-90 -120 -120; 90 130 90];
%     curPos = [1,10,-20];
% %     curPos = [-29,-14,-50]; % 03
% %     curPos = [-3,30,-43]; % 03
% elseif s==3
%     imDims = [-90 -120 -100; 90 130 110];
%     curPos = [1,4,38];
% end

if ~exist('varargin','var') 
    NrFigures = 1;
elseif isempty(varargin)
    NrFigures = 1;
else
    NrFigures = varargin{1};
end

% functionals to ACPC space
% settings:
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = ni.pixdim(1:3);

% functional to ACPC
imgVol = ni.data(:,:,:,1); % do the first
imgVol = imgVol/max(imgVol(:));
[imgSlice1,x1,y1,z1]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
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
    imgSlice=imrotate(imgSlice,90);
elseif sliceThisDim==2
    y1=z1; y=z;
    
    % rotate the images
    imgSlice1=imrotate(imgSlice1,90);
    imgSlice=imrotate(imgSlice,90);
    
elseif sliceThisDim==3
    % x and y stay the same
end

% Make the figure:
if NrFigures > 1
    figure('Position',[0 0 1000 500])
    subplot(1,2,1)
    imagesc(x1,y1,imgSlice1)
    hold on
    axis image
    colormap gray 
else
    figure('Position',[0 0 500 500])
end

if NrFigures>1
    subplot(1,2,2)
end

cm=colormap(jet);

image(x,y,cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:))); % background
hold on
axis image
% tranform imgSlice1 (functionals) to colormap values that I want to use
imgSlice1_color=cat(3,zeros(size(imgSlice1)),zeros(size(imgSlice1)),zeros(size(imgSlice1)));
for k_x = 1:length(x1)
    for k_y = 1:length(y1)
        imgSlice1_color(k_y,k_x,:) = cm(1+floor(abs(imgSlice1(k_y,k_x))*63),:);
    end
end
h=image(x1,y1,imgSlice1_color); % overlay colormap on image
set(h,'AlphaData',.2*ones(size(imgSlice1))) % make transparent

