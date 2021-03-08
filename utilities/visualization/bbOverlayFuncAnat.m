function [imgSlice,x,y,imgSlice1,x1,y1] = bbOverlayFuncAnat(ni,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,varargin)
% Function to plot a functional and overlay with the anatomy

% Inputs: 
%   ni: functional scan
%   niAnatomy: anatomical scan
%   acpcXform: transformation matrix from nifti ijk to anatomy mm
%
%   varargin{1} = maximum value of the overlay.
%
% Outputs:
% imgSlice: anatomical slide
% imgSlice1: overlay in the same space
% 
% There is still an error in case the x, y, z directions differ between the
% anatomy and the functional scans, because the sliceThisDim assumes that
% the directions are the same...
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

if isempty(varargin)
    maxOverlay = max(double(ni.data(:))); % scale to maximum of data
    makeFig = 1;
    minOverlay = [];
else
    if ~isempty(varargin{1})
        maxOverlay = varargin{1}; % scale to maximum of data
    else
        maxOverlay = max(double(ni.data(:))); % scale to maximum of data
    end
    if length(varargin)>1
        if ~isempty(varargin{2})
            makeFig = varargin{2}; % make a figure
        else
            makeFig = 1;
        end
    end
    if length(varargin)>2
        if ~isempty(varargin{3})
            minOverlay = varargin{3}; % scale only show values above 
        else
            minOverlay = [];
        end
    end
end

% slice in ACPC-space to plot
sliceNum = curPos(sliceThisDim);

% Functionals to ACPC space
img2std = acpcXform;
interpType = 'n';
mmPerVox = ni.pixdim(1:3);
imgVol = double(ni.data); 
imgVol = imgVol/maxOverlay; % scale to max
imgVol(imgVol>1) = 1; % set everything larger than max to max
imgVol(imgVol<-1) = -1; % set everything smaller than -max to -max
imgVol = (imgVol+1)/2; % everything is scaled 0:1 now, with zero values at 0.5
[imgSlice1,x1,y1,z1] = dtiGetSlice(img2std,imgVol,sliceThisDim,sliceNum,imDims,interpType,mmPerVox);
if sliceThisDim == 1 || sliceThisDim == 3
    x1 = x1(1,:)';
    y1 = y1(:,1);
elseif sliceThisDim == 2
    x1 = x1(:,1);
    y1 = y1(1,:)';
end
z1 = z1(1,:)';
clear imgVol

% Anatomy to ACPC
img2std = niAnatomy.qto_xyz;
interpType = 'n';
mmPerVox = [1 1 1];
imgVol = niAnatomy.data;
[imgSlice,x,y,z] = dtiGetSlice(img2std,imgVol,sliceThisDim,sliceNum,imDims,interpType,mmPerVox);
if sliceThisDim == 1 || sliceThisDim == 3
    x = x(1,:)';
    y = y(:,1);
elseif sliceThisDim == 2
    x = x(:,1);
    y = y(1,:)';
end
z = z(1,:)';
clear imgVol

% for x and y for plotting:
if sliceThisDim == 1
    x1 = z1; 
    x = z;
    % flip x and y
    x1_t = x1;
    x1 = y1;
    y1 = x1_t;
    x_t = x;
    x = y;
    y = x_t;

    y = y(end:-1:1);
    y1 = y1(end:-1:1);
    
    % and for the images
    imgSlice1 = imrotate(imgSlice1,90);
    imgSlice = imrotate(imgSlice,90);
elseif sliceThisDim == 2
    y1 = z1; 
    y = z;
    
    % rotate the images
    imgSlice1 = imrotate(imgSlice1,90);
    imgSlice = imrotate(imgSlice,90);
    
elseif sliceThisDim == 3
    % x and y stay the same
end

if makeFig==1
    figure('Position',[0 0 500 500])
end

if ~isempty(minOverlay) % remove vals < minOverlay if value
    minTh = [0.5-((minOverlay/maxOverlay)/2) 0.5+((minOverlay/maxOverlay)/2)];
    imgSlice1(imgSlice1>minTh(1) & imgSlice1<minTh(2)) = NaN;
end

subplot(1,5,[1:4])
% cm = colormap(jet);
load loc_colormap
brightness_increase = 3; % factor for increasing brightness of background T1 image
image(x,y,brightness_increase*cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:))); % background
hold on
axis image
% tranform imgSlice1 (functionals) to colormap values that I want to use
imgSlice1_color = cat(3,zeros(size(imgSlice1)),zeros(size(imgSlice1)),zeros(size(imgSlice1)));
for k_x = 1:length(x1)
    for k_y = 1:length(y1)
        if ~isnan(imgSlice1(k_y,k_x)) && (imgSlice1(k_y,k_x)) ~= 0
            imgSlice1_color(k_y,k_x,:) = cm(1+floor(abs(imgSlice1(k_y,k_x))*63),:);
        elseif ~isnan(imgSlice1(k_y,k_x)) && (imgSlice1(k_y,k_x)) == 0
            imgSlice1_color(k_y,k_x,:) = [NaN NaN NaN];
        else
            imgSlice1_color(k_y,k_x,:) = [NaN NaN NaN];
        end
    end
end

h = image(x1,y1,imgSlice1_color); % overlay colormap on image
set(h,'AlphaData',.4*ones(size(imgSlice1))) % make transparent

axis xy

subplot(1,5,5),hold on
cm_vals = 0:maxOverlay/(size(cm,1)-1):maxOverlay;
for kk = 1:length(cm_vals)
    plot([1],cm_vals(kk),'.','Color',cm(kk,:),'MarkerSize',40)
end
text(1.5,cm_vals(end),['>=' num2str(maxOverlay)])
axis off

subplot(1,5,[1:4])