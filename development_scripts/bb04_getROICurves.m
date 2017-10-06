clear all
close all

% Dora Hermes, 2017 


%% The T2* data are here.  
clear all
% close all
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

s_nr = 3;
scan_nr = 3;
s_info = bb_subs(s_nr);
subj=s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii']));

% % Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));
% % load Veno rotation matrix:
% xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));

% now do an SVD on the odd responses:

% load PPG responses
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};

% nifti:
ni=niftiRead(fullfile(dDir,subj,scan,[scanName '.nii.gz']));
% load coregistration matrix (for the functionals):
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new;

load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT']),'t')

% load average of all odd heartbeats:
ppgTS=niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_odd.nii.gz']));

% load average of all odd heartbeats:
ppgTSeven=niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_even.nii.gz']));

% Load the correlation with heartbeat (made with bbCorrelate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

% Scale the time-series matrix by the COD
% divide by max
ppgTS.data = ppgTS.data ./ repmat(max(abs(ppgTS.data),[],4),[1,1,1,size(ppgTS.data,4)]);
ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
% multiply by correlation size (absolute)
ppgTS.data = ppgTS.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTS.data,4)]));
ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));

% reshape for SVD
a = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);

a = a(:,t>=-.5 & t<=2);
t_sel = t(t>=-.5 & t<=2);
a = a-repmat(mean(a,2),1,size(a,2)); % subtract the mean
[u,s,v]=svd(a','econ');
s=diag(s);

%%%% test for the sign of the 2nd component (also check for 1st???)
% in the 2nd component, the first peak should be negative
[~,pm_i]=findpeaks(double(u(:,2)),'minpeakdistance',10);
[~,nm_i]=findpeaks(-double(u(:,2)),'minpeakdistance',10);
first_peak = min([pm_i; nm_i]);
if u(first_peak,2)<0
    % do nothing
elseif u(first_peak,2)>0
    % reverse sign pc and weights
    u(:,2)=-u(:,2); v(:,2) = -v(:,2);
end

%%%% put output in structures:
% put 2 components weights in a matrix
out = [];
for k=1:2
    out(k).weights = reshape(v(:,k),[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);
end

% %%%%% MODEL WITH 2 COMPONENTS:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';
svdResults.model = reshape(pred,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3) size(ppgTS.data,4)]);

%%
%% Get ROI indices in the functional scan 
%%

% ROInames = {'R_lat_ventr','R_inferior_lat_ventr','R_choroid_plexus','3rd_ventr','4th_ventr','CSF'};
ROInames = {'R_choroid_plexus','R_lat_ventr','R_inferior_lat_ventr','L_choroid_plexus','L_lat_ventr','L_inferior_lat_ventr','CSF','3rd_ventr','4th_ventr'};
ROIplotInd = [1 3 5 2 4 6 7 9 11];
figure('Position',[0 0 450 700])

for rr = 1:length(ROInames)

    niROIname = ROInames{rr};
    niROI = niftiRead(fullfile(dDir,subj,'freesurfer','nii',[niROIname '.nii.gz']));

    % get ROI indices:
    [xx,yy,zz] = ind2sub(size(niROI.data),find(niROI.data>0));

    % now ROI indices to ACPC (mm):
    xyz_acpc = mrAnatXformCoords(niROI.qto_xyz, [xx,yy,zz]);
    clear xx yy zz % housekeeping

    % now ACPC coordinates to functional indices:
    ijk_func = mrAnatXformCoords(inv(acpcXform), xyz_acpc);
    ijk_func = round(ijk_func); % round to get indices
    ijk_func = unique(ijk_func,'rows'); % only take unique voxels

%     %%%% check for coordinates in functional space
%     z_slice = [5 10 15 18 20 23 25 26 27 30 33 36];
%     figure
%     for kk = 1:length(z_slice)       
%         subplot(3,4,kk)
%         imagesc(ni.data(:,:,z_slice(kk))),hold on
%         xyz_plot = ijk_func(ijk_func(:,3)==z_slice(kk),:);
%         plot(xyz_plot(:,2),xyz_plot(:,1),'r.')
%         axis image 
%     end
    % get ROI curves
    % imgVol = out(1).weights; % PC1 weight
    % imgVol = out(2).weights; % PC2 weight
%     imgVol = svdResults.model;
    imgVol = ppgTS.data;

    roiTrace = NaN(size(ijk_func,1),size(imgVol,4));

    for kk = 1:size(ijk_func,1)
        roiTrace(kk,:) = imgVol(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3),:);
        % note that x and y are only switched around for plotting position, 
        % not for getting the actual image values
    end

    subplot(length(ROInames)-3,2,ROIplotInd(rr)),hold on
    
    plot(t,roiTrace')
%     trace2plot = mean(roiTrace,1);
%     trace2plot_std = 2*std(roiTrace,[],1)/sqrt(size(roiTrace,1));

%     trace2plot = median(roiTrace,1);
%     trace2plotPlusErr = quantile(roiTrace,.84,1);
%     trace2plotMinErr = quantile(roiTrace,.16,1);   
%     fill([t t(end:-1:1)],[trace2plotMinErr trace2plotPlusErr(end:-1:1)],[.7 .7 .7],'EdgeColor',[.7 .7 .7])
%     plot(t,trace2plot,'k','LineWidth',2)
%     plot([0 0],[min(trace2plotMinErr) max(trace2plotPlusErr)],'k')
    
    ylabel(niROIname)
    xlim([t(1) 1.5])
    
    subplot(length(ROInames)-3,2,2*(length(ROInames)-3))
    title('mean model and 68% quantiles')
end

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/ROI/' subj '_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_ROIsTest'])



%%
%% Plot one ROI at the time 
%%

% ROInames = {'R_choroid_plexus','R_lat_ventr','R_inferior_lat_ventr','L_choroid_plexus','L_lat_ventr','L_inferior_lat_ventr','CSF','3rd_ventr','4th_ventr'};
ROInames = {'R_lat_ventr'};

for rr = 1:length(ROInames)

    niROIname = ROInames{rr};
    niROI = niftiRead(fullfile(dDir,subj,'freesurfer','nii',[niROIname '.nii.gz']));

    % get ROI indices:
    [xx,yy,zz] = ind2sub(size(niROI.data),find(niROI.data>0));

    % now ROI indices to ACPC (mm):
    xyz_acpc = mrAnatXformCoords(niROI.qto_xyz, [xx,yy,zz]);
    clear xx yy zz % housekeeping

    % now ACPC coordinates to functional indices:
    ijk_func = mrAnatXformCoords(inv(acpcXform), xyz_acpc);
    ijk_func = round(ijk_func); % round to get indices
    ijk_func = unique(ijk_func,'rows'); % only take unique voxels

%     %%%% check for coordinates in functional space
%     z_slice = [5 10 15 18 20 23 25 26 27 30 33 36];
%     figure
%     for kk = 1:length(z_slice)       
%         subplot(3,4,kk)
%         imagesc(ni.data(:,:,z_slice(kk))),hold on
%         xyz_plot = ijk_func(ijk_func(:,3)==z_slice(kk),:);
%         plot(xyz_plot(:,2),xyz_plot(:,1),'r.')
%         axis image 
%     end
    % get ROI curves
    % imgVol = out(1).weights; % PC1 weight
    % imgVol = out(2).weights; % PC2 weight
%     imgVol = svdResults.model;
    imgVol = ppgTS.data;
    % Mask by reliability R
%     imgVol(ppgR.data<.5) = NaN;
    
    roiTrace = NaN(size(ijk_func,1),size(imgVol,4));
    svdVals = NaN(size(ijk_func,1),2);

    for kk = 1:size(ijk_func,1)
        roiTrace(kk,:) = imgVol(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3),:);
        % note that x and y are only switched around for plotting position, 
        % not for getting the actual image values
        svdVals(kk,:) = [out(1).weights(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3)) out(2).weights(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3))];
    end
    % Remove voxels with low R
    roiTrace(isnan(roiTrace(:,1)),:) = [];

    figure('Position',[0 0 800 500])
    subplot(1,2,1)
    %%%% Imagesc:
    imagesc(t,[1:size(roiTrace,1)],roiTrace,[-1 1])
    
    %%%% Plot all traces:
%     hold on
%     plot(t,roiTrace')
%     plot([0 0],[min(roiTrace(:)) max(roiTrace(:))],'k')

    %%%% Mesh or surf
%     xx = repmat([1:size(roiTrace,1)],length(t),1)';
%     yy = repmat(t,size(roiTrace,1),1);
%     mesh(xx,yy,roiTrace)
%     view(100,30)

%     xlim([t(1) 1.5])
    
    title(['Traces sub ' int2str(s_nr) '  scan ' int2str(scan_nr) ' roi: ' niROIname])
    subplot(1,2,2)
    
    title('SVD PC1 and PC2 weights')
    imagesc(svdVals,[-.02 .02])
end

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/ROI/' subj '_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_ROIsTest'])

% Get coordinates back in acpc:
xyz_acpc_sparse = mrAnatXformCoords(acpcXform, ijk_func);