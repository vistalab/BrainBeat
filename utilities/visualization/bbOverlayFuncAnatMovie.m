function [] = bbOverlayFuncAnatMovie(ni,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,videoName,t)
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

% functionals to ACPC space
% settings:
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = ni.pixdim(1:3);

% functional to ACPC
imgVol = ni.data; % do the first
% imgVol = imgVol/max(imgVol(:));
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


%%
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

%%
fid = figure('Position',[0 0 500 500]);

vidObj = VideoWriter(videoName,'MPEG-4'); %

open(vidObj); 

for kk = 1:size(ni.data,4) % loop over 4th dimension in the nifti
    % Background:
    imagesc(x,y,imgSlice/max(imgSlice(:)));
    set(gca,'CLim',[0 .3])
    hold on
    axis image
    colormap gray
        
    % Overlay:
    currentSlice = imgSlice1(:,:,kk);
    % Initialize imgSlice1 and transform to colormap values that I want to use
    imgSlice1_color=cat(3,zeros(size(currentSlice)),zeros(size(currentSlice)),zeros(size(currentSlice)));
    for k_x = 1:length(x1)
        for k_y = 1:length(y1)
            imgSlice1_color(k_y,k_x,:) = cm(floor(currentSlice(k_y,k_x)*99)+100,:);
        end
    end
    h = imagesc(x1,y1,imgSlice1_color,[0 1]); % overlay colormap on image
    set(h,'AlphaData',.2*ones(size(imgSlice1,1),size(imgSlice1,2))) % make transparent
    
    title(['t = ' num2str(t(kk),3)])
    
    % Write each frame to the file.
    for m=1:5 % write X frames: decides speed
        writeVideo(vidObj,getframe(fid));
    end
    
    clf
    
end

close(vidObj);

