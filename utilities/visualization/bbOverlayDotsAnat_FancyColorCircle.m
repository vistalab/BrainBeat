function [] = bbOverlayDotsAnat_FancyColorCircle(ni1,ni2,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,varargin)
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
%
% Dora Hermes, 2017 

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
    maxPlot = max(ni1.data(:));
    pointSize = 12;
    dotTh = NaN;
else
    if length(varargin)==1
        maxPlot = varargin{1};
        pointSize = 12;
        dotTh = NaN;
    elseif length(varargin)==2
        maxPlot = varargin{1};
        pointSize = varargin{2};
        dotTh = NaN;
    elseif length(varargin)==3
        maxPlot = varargin{1};
        pointSize = varargin{2};
        dotTh = varargin{3};
    end
end

% functionals to ACPC space
% settings:
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = ni1.pixdim(1:3);

% functional to ACPC
imgVol1 = ni1.data;
imgVol1 = imgVol1/maxPlot; % rescale: and set max from -1:1
imgVol1(imgVol1>1) = 1; % clip larger vals
imgVol1(imgVol1<-1) = -1; % clip smaller vals
imgVol2 = ni2.data;
imgVol2 = imgVol2/maxPlot;% rescale: and set max from -1:1
imgVol2(imgVol2>1) = 1; % clip larger vals
imgVol2(imgVol2<-1) = -1; % clip smaller vals
% get the slice from the Color
[imgSlice1,x1,y1,z1] = dtiGetSlice(img2std, imgVol1, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
% get the slice from the Intensity
[imgSlice2] = dtiGetSlice(img2std, imgVol2, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
if sliceThisDim == 1 || sliceThisDim == 3
    x1 = x1(1,:)';
    y1 = y1(:,1);
elseif sliceThisDim == 2
    x1 = x1(:,1);
    y1 = y1(1,:)';
end
z1 = z1(1,:)';

% Anatomy to ACPC
imgVol = niAnatomy.data;
img2std = niAnatomy.qto_xyz;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [1 1 1];
[imgSlice,x,y,z] = dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
if sliceThisDim == 1 || sliceThisDim == 3
    x = x(1,:)';
    y = y(:,1);
elseif sliceThisDim == 2
    x = x(:,1);
    y = y(1,:)';
end
z = z(1,:)';

% for x and y for plotting:
if sliceThisDim==1
    x1 = z1; x = z;
    % flip x and y
    x1_t = x1;
    x1 = y1;
    y1 = x1_t;
    x_t = x;
    x = y;
    y = x_t;
    
    % and for the images
    imgSlice1 = imrotate(imgSlice1,90);
    imgSlice2 = imrotate(imgSlice2,90);
    imgSlice = imrotate(imgSlice,90);
elseif sliceThisDim==2
    y1 = z1; y = z;
    
    % rotate the images
    imgSlice1 = imrotate(imgSlice1,90);
    imgSlice2 = imrotate(imgSlice2,90);
    imgSlice = imrotate(imgSlice,90);
    
elseif sliceThisDim==3
    % x and y stay the same
end

% now plot stuff:
figure('Position',[0 0 600 500])

colormap gray 

subplot(1,1,1)
% Background:
imagesc(x,y,imgSlice/max(imgSlice(:)));
colormap gray
set(gca,'CLim',[0 .3])
hold on
axis image
 
% Dots on image:
for k = 1:size(imgSlice1,1)
    for m = 1:size(imgSlice1,2)
        % get plotting location
        k_x = y1(k);
        m_y = x1(m);

        v1 = imgSlice1(k,m);
        v2 = imgSlice2(k,m);
        
        if isnan(dotTh) % just plot everything
            data_colors_rgb = bbData2Colors([v1,v2]);
            plot(m_y+1,k_x+1,'.','Color',data_colors_rgb,'MarkerSize',pointSize)
        else
            if abs(v1)>dotTh || abs(v2)>dotTh
                data_colors_rgb = bbData2Colors([v1,v2]);
                plot(m_y+1,k_x+1,'.','Color',data_colors_rgb,'MarkerSize',pointSize)
            end
        end
    end
end
%% to plot colorscale
% 
% % make 2D colormap: one to vary color, the other varies intensity
% cm = jet(250); cm = cm(26:225,:);
% cm = cm(end:-1:1,:);
% cm = cm+.4; cm(cm>1)=1;
% gray_vect = .2*ones(200,3);
% cm2D = zeros(100,size(cm,1),3);
% for kk = 1:100
%     cm2D(kk,:,:) = cm*kk/100 + gray_vect*(100-kk)/100;
% end
% 
% figure,hold on
% for kk = 1:size(cm2D,1)
%     for mm = 1:size(cm2D,2)
%         plot(kk,mm,'.','MarkerSize',25,'Color',cm2D(kk,mm,:))
%     end
% end
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-depsc',[dDir './figures/2DcolormapNew'])
% print('-painters','-r300','-dpng',[dDir './figures/2DcolormapNew'])

