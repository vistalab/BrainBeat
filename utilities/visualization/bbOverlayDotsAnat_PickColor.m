function [] = bbOverlayDotsAnat_PickColor(ni,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,cm,maxPlot)
% function to plot a functional and overlay with the anatomy

% Inputs: 
%   ni: functional scan
%   niAnatomy: anatomical scan
%   acpcXform: transformation matrix from nifti ijk to anatomy mm
% 
% Examples:
%
% This example scales the dots to the maximum from the image:
%   bbOverlayDotsAnat(niOverlay,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,cm)
%
% %%% now overlay r-map with anatomy

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

% functionals to ACPC space
% settings:
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = ni.pixdim(1:3);

% functional to ACPC
imgVol = ni.data;

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
% figure('Position',[0 0 600 500])
% cm=hsv(200);

colormap gray 

subplot(1,5,1:4)

% Plot background:
imagesc(x,y,imgSlice/max(imgSlice(:)));
colormap gray
set(gca,'CLim',[0 1])
hold on
axis image

imgSlice1 = round(length(cm)*(imgSlice1./maxPlot));

% Plot dots on image
for k = 1:size(imgSlice1,1)
    for m = 1:size(imgSlice1,2)
        % get plotting location
        k_x = y1(k);
        m_y = x1(m);

        val_plot = imgSlice1(k,m);
        
        if val_plot>0
            c_use = cm(val_plot,:); % fot plotting plus and min
            plot(m_y+1,k_x+1,'.','Color',c_use,'MarkerSize',8)
        end        
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
drawnow

subplot(1,5,1:4)
