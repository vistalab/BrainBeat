function [voxelTimeseries] = bbGetVoxelTimeseries(ni,acpcXform,curPos,varargin)
% function to plot a functional and overlay with the anatomy

% Inputs: 
%   ni: functional scan, plots the mean if a series
%   acpcXform: transformation matrix from nifti ijk to anatomy mm
%   curpos: voxel to get timeseries from in mm 

xyz_acpc = curPos; %ACPC coordinates

% ACPC coordinates to functional indices:
ijk_func = mrAnatXformCoords(inv(acpcXform), xyz_acpc);
ijk_func = round(ijk_func); % round to get indices
ijk_func = unique(ijk_func,'rows'); % only take unique voxels

%     %%%% check for coordinates in functional space
%     z_slice = ijk_func(3);
%     figure
%     imagesc(mean(ni.data(:,:,z_slice,1),4)),hold on
%     xyz_plot = ijk_func(ijk_func(:,3)==z_slice,:);
%     plot(xyz_plot(:,2),xyz_plot(:,1),'r.') % note that x/y is only
%     shifted for plotting
%     axis image 
    % get ROI curves
    % imgVol = out(1).weights; % PC1 weight
    % imgVol = out(2).weights; % PC2 weight
%     imgVol = svdResults.model;
imgVol = ni.data;

voxelTimeseries = squeeze(imgVol(ijk_func(1),ijk_func(2),ijk_func(3),:));
% note here that only an xy /switch for plotting is necessary, but not to get the actual curves!!!


%% This code also works, but is much longer ...

% if ~isempty(varargin)
%     ifsliceThisDim = varargin{1};
%     imDims = varargin{2};
% end

% % functionals to ACPC space
% % settings:
% img2std = acpcXform;
% sliceNum = curPos(sliceThisDim);
% interpType = 'n';
% mmPerVox = ni.pixdim(1:3);
% 
% % functional to ACPC
% imgVol = ni.data; % plot the mean timeseries
% [imgSlice1,x1,y1,z1]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
% if sliceThisDim == 1 || sliceThisDim == 3
%     x1=x1(1,:)';
%     y1=y1(:,1);
% elseif sliceThisDim == 2
%     x1=x1(:,1);
%     y1=y1(1,:)';
% end
% z1=z1(1,:)';
% 
% % for x and y for plotting:
% if sliceThisDim==1
%     x1=z1; x=z;
%     % flip x and y
%     x1_t=x1;
%     x1=y1;
%     y1=x1_t;
%     x_t=x;
%     x=y;
%     y=x_t;
%     
%     % and for the images
%     imgSlice1=imrotate(imgSlice1,90);
% elseif sliceThisDim==2
%     y1=z1; y=z;
%     
%     % rotate the images
%     imgSlice1=imrotate(imgSlice1,90);
%     
% elseif sliceThisDim==3
%     % x and y stay the same
% end
% 
% %% Get the timeseries for a voxel 
% 
% [~,x_ind] = min(abs(curPos(1)-y1));
% [~,y_ind] = min(abs(curPos(2)-x1));
%          
% voxelTimeseries = squeeze(imgSlice1(x_ind,y_ind,:));

