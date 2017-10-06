function [] = bbOverlayTimeseriesVeno(ni,niColor,niAnatomy,acpcXform,acpcXformVeno,sliceThisDim,imDims,curPos)
% function to plot a functional and overlay with the anatomy

% Inputs: 
%   ni: functional scan, plots the mean if a series
%   ni_time: sample time of the timeseries in the nifti, used for scaling
%   niAnatomy: anatomical scan
%   acpcXform: transformation matrix from nifti ijk to anatomy mm
%
%   The plotted X and Y do not correspond to actual xyz/mm coordinates, but
%   are in a different frame. The input XYZ are in mm coordinates

% functionals to ACPC space
% settings:
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = ni.pixdim(1:3);

% functional to ACPC
imgVol = ni.data; % plot the mean timeseries
[imgSlice1,x1,y1,z1]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
% colors to ACPC
imgVolColor = niColor.data; % plot the mean timeseries
[imgSlice1Color]=dtiGetSlice(img2std, imgVolColor, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
if sliceThisDim == 1 || sliceThisDim == 3
    x1=x1(1,:)';
    y1=y1(:,1);
elseif sliceThisDim == 2
    x1=x1(:,1);
    y1=y1(1,:)';
end
z1=z1(1,:)';

% Veno to ACPC
imgVol = niAnatomy.data;
img2std = acpcXformVeno;
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
    imgSlice1Color=imrotate(imgSlice1Color,90);
    imgSlice=imrotate(imgSlice,90);
elseif sliceThisDim==2
    y1=z1; y=z;
    
    % rotate the images
    imgSlice1=imrotate(imgSlice1,90);
    imgSlice1Color=imrotate(imgSlice1Color,90);
    imgSlice=imrotate(imgSlice,90);
    
elseif sliceThisDim==3
    % x and y stay the same
end

%%
%%%%% overlay with time series
figure('Position',[0 0 1000 600])
cm=colormap(jet);
subplot(1,2,1)
% show the background:
imagesc(x1,y1,imgSlice1(:,:,1),[0 1]);
hold on
colormap gray
axis image
title('First time point in nifti')

subplot(1,2,2)
% show the background:
%     image(x,y,cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:)));
imagesc(x,y,imgSlice/max(imgSlice(:)),[0 1]);
% set(gca,'CLim',[0 .3])
hold on
axis image

for sub_p=1:2
    subplot(1,2,sub_p)
    % Get scale factor the curves:
    max_plot = max(abs(ni.data(:)));
    
    % Get the curve colors:
    imgColors = imgSlice1Color; % R^2 has a maximum of 1, this should be maintained
    
    for k=1:size(imgSlice1,1) % rows
        for m=1:size(imgSlice1,2) % columns
            % set the max of the curves to 1
            r_curve=squeeze(imgSlice1(k,m,:))./max_plot;

            % get plotting location: index k,m is plotted at m,k
            k_x = y1(k);
            m_y = x1(m);

            % scaling factor to make the curve fit in a 4 mm voxels 
            s_f=3; 

            % select the color for the curve:
            c_use=cm(1+floor(imgColors(k,m)*63),:);
            
            % This makes sure that the curves have the image top as up,
            % even if the y-axis is reversed
            if isequal(get(gca,'YDir'),'reverse')
                plot(m_y+s_f*[0:size(imgSlice1,3)-1]/size(imgSlice1,3)-2,...
                    k_x-s_f*r_curve,...
                    'Color',c_use)
            else
                plot(m_y+s_f*[0:size(imgSlice1,3)-1]/size(imgSlice1,3)-2,...
                    k_x+s_f*r_curve,...
                    'Color',c_use)
            end
    %         plot(m_y-2,k_x,'.','Color',[.5 .5 .5])
        end
    end
end


