function [] = bbOverlayDotsAnatMovie(ni,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,videoName,t)
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

cm1 = spring(90);
a = [(0:1/9:1)' zeros(10,1) (0:1/9:1)']; 
cm1 = [a; cm1];
cm2 = winter(100);
a = [zeros(10,1) zeros(10,1) (0:1/9:1)']; 
cm2 = [a; cm2];
cm = [cm1(end:-1:1,:); cm2];

% cm = hot(210);

%% Upsample timeseries (imgSlice1)

upFactor = 4;
tUp = interp(t,upFactor);
imgSliceUp = zeros(size(imgSlice1,1),size(imgSlice1,2),length(tUp));

for kx = 1:size(imgSlice1,1)
    for ky = 1:size(imgSlice1,2)
        imgSliceUp(kx,ky,:) = interp(squeeze(imgSlice1(kx,ky,:)),upFactor);
        %imgSliceUp(kx,ky,:) = smooth(imgSliceUp(kx,ky,:),20); %maybe smooth
    end
end
    

%%
fid = figure('Position',[0 0 500 500]);

vidObj = VideoWriter(videoName,'MPEG-4'); %

open(vidObj); 

for kk = 1:size(imgSliceUp,3) % loop over time
    % Background:
    imagesc(x,y,imgSlice/max(imgSlice(:)));
    set(gca,'CLim',[0 .3])
    hold on
    axis image
    colormap gray
    
    % Overlay:
    currentSlice = imgSliceUp(:,:,kk);
    % Plot dots on image
    for kx = 1:size(currentSlice,1)
        for mx = 1:size(currentSlice,2)
            % get plotting location
            k_x = y1(kx);
            m_y = x1(mx);

            val_plot = currentSlice(kx,mx);
            if val_plot > 1
                val_plot = 1;
            elseif val_plot < -1
                val_plot = -1;
            end

            c_use = cm(ceil(val_plot*99.5)+100,:); % fot plotting plus and min

            plot(m_y+1,k_x+1,'.','Color',c_use,'MarkerSize',8)
        end
    end
    
    title(['t = ' num2str(tUp(kk),3)])
    
    % Let the number of frames to write depend on the timing, such that
    % peaks can be emphasized
    nr_frames = 1;
%     if t(kk)>-0.3 && t(kk)<0.3 
%         nr_frames = 20;
%     else
%         nr_frames = 5;
%     end
    
    % Write each frame to the file
    for m = 1:nr_frames % write X frames: decides speed
        writeVideo(vidObj,getframe(fid));
    end
    
    clf
    
end

close(vidObj);

