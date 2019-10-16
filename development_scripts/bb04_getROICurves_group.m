clear all
close all

% Dora Hermes, 2017 

%% The T2* data are here.  
clear all
% close all
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

s_nr = 4;
scan_nr = 3;
s_info = bb_subs(s_nr);
subj=s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

% % Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));
% % load Veno rotation matrix:
% xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));

% now do an SVD on the odd responses:

% load PPG responses
scan = s_info.scan{scan_nr};
scanName = s_info.scanName{scan_nr};

% nifti:
ni = niftiRead(fullfile(dDir,subj,scan,[scanName '.nii.gz']));
% load coregistration matrix (for the functionals):
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new;

load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT']),'t')

% load all odd heartbeats:
ppgTS = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse.nii.gz']));

% load average of all odd heartbeats:
ppgTSodd = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_odd.nii.gz']));

% load average of all odd heartbeats:
ppgTSeven = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_even.nii.gz']));

% Load the correlation with heartbeat (made with bbCorrelate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

% Scale the time-series matrix by the COD
% divide by max
ppgTSodd.data = ppgTSodd.data ./ repmat(max(abs(ppgTSodd.data),[],4),[1,1,1,size(ppgTSodd.data,4)]);
ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
% multiply by correlation size (absolute)
ppgTSodd.data = ppgTSodd.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSodd.data,4)]));
ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));

% reshape for SVD
a = reshape(ppgTSodd.data,[numel(ppgTSodd.data(:,:,:,1)) length(t)]);

a = a(:,t>=-.5 & t<=2);
t_sel = t(t>=-.5 & t<=2);
a = a-repmat(mean(a,2),1,size(a,2)); % subtract the mean
[u,s,v] = svd(a','econ');
s=diag(s);

%%%% test for the sign of the 2nd component (also check for 1st???)
% in the 2nd component, the first peak should be negative
[~,pm_i] = findpeaks(double(u(:,2)),'minpeakdistance',10);
[~,nm_i] = findpeaks(-double(u(:,2)),'minpeakdistance',10);
first_peak = min([pm_i; nm_i]);
if u(first_peak,2)<0
    % do nothing
elseif u(first_peak,2)>0
    % reverse sign pc and weights
    u(:,2) = -u(:,2); v(:,2) = -v(:,2);
end

%%%% put output in structures:
% put 2 components weights in a matrix
out = [];
for k = 1:2
    out(k).weights = reshape(v(:,k),[size(ppgTSodd.data,1) size(ppgTSodd.data,2) size(ppgTSodd.data,3)]);
end

% %%%%% MODEL WITH 2 COMPONENTS:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';
svdResults.model = reshape(pred,[size(ppgTSodd.data,1) size(ppgTSodd.data,2) size(ppgTSodd.data,3) size(ppgTSodd.data,4)]);

%%
%% Get ROI indices in the functional scan 
%%

% load freesurfer segmentation (these are resliced to functional space)
niFs = niftiRead(fullfile(dDir,subj,scan,[scanName '_r_aseg_auto.nii.gz']));
% FS_rois = [...
%         3 % left gray matter
%         42% right gray matter
%         2 % left white matter
%         41 % right white matter
%         4 % left lateral ventricle
%         5 % left inferior lateral ventricle
%         14 % 3rd ventricle
%         15 % 4th ventricle
%         24 % CSF
%         31 % left choroid plexus
%         43 % right lateral ventricle
%         44 % right inferior lateral ventricle
%         63 % right choroid plexus
%         72]; % 5th ventricle

% get manually segmented ROIs - these are in anatomical space
clear niROI
roi_list = {'CFlowvoids','AnteriorSSS','SSS','LeftTransverse','RightTransverse'};
for rr = 1:length(roi_list)
    niROIname       = roi_list{rr};
    niROI(rr).ni    = niftiRead(fullfile(dDir,subj,s_info.anat,['r' s_info.anatName niROIname '.nii']));
end

% create output structure with timeseries per ROI
clear roiTS

% get timeseries for all manually segmented ROIs
for rr = 1:length(roi_list)

    % get ROI indices:
    [xx,yy,zz] = ind2sub(size(niROI(rr).ni.data),find(niROI(rr).ni.data>.5));

    % now ROI indices to ACPC (mm):
    xyz_acpc = mrAnatXformCoords(niROI(rr).ni.qto_xyz, [xx,yy,zz]);
    clear xx yy zz % housekeeping

    % now ACPC coordinates to functional indices:
    ijk_func = mrAnatXformCoords(inv(acpcXform), xyz_acpc);
    ijk_func = round(ijk_func); % round to get indices
    ijk_func = unique(ijk_func,'rows'); % only take unique voxels
%     ijk_func(min(ijk_func,[],2)<0,:) = []; % empty voxels outside of functional

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
    
    imgVol1 = 100*ppgTS.data; %timeseries in percent signal change
    imgVol2 = svdResults.model;
    
    roiTS(rr).roiCod = NaN(size(ijk_func,1),1);
    roiTS(rr).roiWeights = NaN(size(ijk_func,1),2);
    roiTS(rr).roiTrace = NaN(size(ijk_func,1),size(imgVol1,4));
    roiTS(rr).pred = NaN(size(ijk_func,1),size(imgVol2,4));

    for kk = 1:size(ijk_func,1)
        % note that x and y are only switched around for plotting position, 
        % not for getting the actual image values
        roiTS(rr).roiTrace(kk,:) = imgVol1(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3),:);
        roiTS(rr).pred(kk,:) = imgVol2(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3),:);       
        roiTS(rr).roiCod(kk) = ppgR.data(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3));
        roiTS(rr).roiWeights(kk,1) = out(1).weights(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3));
        roiTS(rr).roiWeights(kk,2) = out(2).weights(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3));

    end
end

%%%% Add freesurfer ROI traces (e.g. CSF)
fs_segm = {[3 42],[2 41],[31 63],[14],[15],[24]};
roi_list = {'CFlowvoids','AnteriorSSS','SSS','LeftTransverse','RightTransverse',...
    'Gray','White','ChoroidPlexus','3rdVentr','4thVentr','CSF'};

% reshape timeseries for voxel selection
tsVect1 = 100*reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);
% reshape prediction for voxel selection
tsVect2 = reshape(svdResults.model,[numel(svdResults.model(:,:,:,1)) length(t)]);

% which indices are added in the output: roiTS
out_ind = 6:5+length(fs_segm);
for rr = 1:length(out_ind)
    % get the voxels from the current ROI:
    thisRoi_voxels = find(ismember(niFs.data,fs_segm{rr}));
    roiTS(out_ind(rr)).Cod = ppgR.data(thisRoi_voxels);
    roiTS(out_ind(rr)).roiWeights(:,1) = out(1).weights(thisRoi_voxels);
    roiTS(out_ind(rr)).roiWeights(:,2) = out(2).weights(thisRoi_voxels);    
    roiTS(out_ind(rr)).roiTrace = tsVect1(thisRoi_voxels,:);
    roiTS(out_ind(rr)).pred = tsVect2(thisRoi_voxels,:);
end
clear tsVect1 tsVect2

%% quick figure
figure('Position',[0 0 300 800])

for kk = 1:length(roiTS)
    subplot(length(roiTS),2,kk*2-1)
    plot(t,mean(roiTS(kk).roiTrace,1),'k')
    subplot(length(roiTS),2,kk*2)
    plot(t,mean(roiTS(kk).pred,1),'r')
end

%% save subject timeseries

tsSaveName = fullfile(dDir,subj,scan,[scanName '_ROItimeseries01.mat']);

save(tsSaveName,'roiTS','t','roi_list')

%%

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/ROI/TimeSerie_sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_roi-CSF_COD-0_3'])
% print('-painters','-r300','-depsc',[dDir './figures/ROI/TimeSerie_sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_roi-CSF_COD-0_3'])

