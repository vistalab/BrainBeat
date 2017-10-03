function [] = bbOverlayDotsAnat(ni,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,varargin)
% function to plot a functional and overlay with the anatomy

% Inputs: 
%   ni: functional scan
%   niAnatomy: anatomical scan
%   acpcXform: transformation matrix from nifti ijk to anatomy mm
% 
% Examples:
%
% This example scales the dots to the maximum from the image:
%   bbOverlayDotsAnat(niOverlay,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
%
% This example scales the dots to the maximum specifed in maxPlot:
%   maxPlot = .1;
%   bbOverlayDotsAnat(niOverlay,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,maxPlot)

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
    maxPlot = max(ni.data(:));
else
    maxPlot = varargin{1};
end

% functionals to ACPC space
% settings:
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = ni.pixdim(1:3);

% functional to ACPC
imgVol = ni.data;
% rescale: and set max from -1:1
imgVol = imgVol/maxPlot;
imgVol(imgVol>1) = 1;
imgVol(imgVol<-1) = -1;
% get the slice
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

% now plot stuff:
figure('Position',[0 0 600 500])
% cm=hsv(200);

% make my colormap: cyan - blue - black - red - yellow
cm1(1:100,2)=[0:1/99:1]';
cm1=[repmat([0 0 0],100,1)];
cm1(1:10,1)=[0:0.7/9:0.7]';
cm1(10:100,1)=[0.7:(1-0.7)/90:1]';
cm1(1:10,2)=[0]';
cm1(1:10,3)=[0]';
cm1(20:100,2)=[0:1/80:1]';
cm2=[repmat([0 0 0],100,1)];
cm2(1:10,3)=[0:0.7/9:0.7]';
cm2(10:100,3)=[0.7:(1-0.7)/90:1]';
cm2(1:10,2)=[0]';
cm2(1:10,1)=[0]';
cm2(20:100,2)=[0:1/80:1]';
cm2=cm2(end:-1:1,:);
cm=[cm2; cm1];

colormap gray 

subplot(1,5,1:4)
image(x,y,cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:))); % background
hold on
axis image
% % tranform imgSlice1 (functionals) to colormap values that I want to use
% imgSlice1_color=cat(3,zeros(size(imgSlice1)),zeros(size(imgSlice1)),zeros(size(imgSlice1)));
% for k_x = 1:length(x1)
%     for k_y = 1:length(y1)
%         imgSlice1_color(k_y,k_x,:) = cm(1+floor(abs(imgSlice1(k_y,k_x))*63),:);
%     end
% end
% h=image(x1,y1,imgSlice1_color); % overlay colormap on image
% set(h,'AlphaData',.2*ones(size(imgSlice1))) % make transparent 

for k=1:size(imgSlice1,1)
    for m=1:size(imgSlice1,2)
        % get plotting location
        k_x = y1(k);
        m_y = x1(m);

        val_plot = imgSlice1(k,m);
        
        c_use=cm(ceil(val_plot*99.5)+100,:); % fot plotting plus and min

        plot(m_y+1,k_x+1,'.','Color',c_use,'MarkerSize',8)
    end
end
subplot(1,5,5),hold on
for kk = 1:size(cm,1)
    plot(1,kk,'.','Color',cm(kk,:),'MarkerSize',20)
end
text(1.5,size(cm,1),'max')
text(1.5,size(cm,1)/2,'0')
text(1.5,0,'min')
axis off

