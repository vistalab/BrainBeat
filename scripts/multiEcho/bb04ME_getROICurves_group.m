clear all
close all

% Dora Hermes, 2017 

%% The T2* data are here.  
% clear all
% close all
% dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% group with FA of 48 includes sub/scan: [1/3, 2/3, 3/3, 4/1, 5/1, 6/1]
% group with FA of 48 includes second sub/scan: [6/2]

s_nr = 4; 
s_info = bb_subs(s_nr);
subj = s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

% get multi echo data
if s_nr==4
    scan_nrs = {[4 5],[6 7]};
elseif s_nr==5
    scan_nrs = {[4 5],[6 7]};
elseif s_nr==6
    scan_nrs = {[4 5],[6 7],[9 10]};
end

scn_me = 1;%:length(scan_nrs)
% first echo
scan1 = s_info.scan{scan_nrs{scn_me}(1)};
scanName1 = s_info.scanName{scan_nrs{scn_me}(1)};

% load coregistration matrix:
load(fullfile(dDir,subj,scan1,[scanName1 'AcpcXform_new.mat']));
acpcXform = acpcXform_new;

% load time T:
t_me = load(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponseT']),'t');

% load average of all heartbeats:
fname_s0 = [scanName1 '_PPGtrigResponse_S0.nii.gz'];
fname_t2s = [scanName1 '_PPGtrigResponse_T2s.nii.gz'];
ni_s0 = niftiRead(fullfile(dDir,subj,scan1,fname_s0));
ni_t2s = niftiRead(fullfile(dDir,subj,scan1,fname_t2s));

% Load the correlation with heartbeat (made with bbCorrelate2physio):
ppgRname_S0 = fullfile(dDir,subj,scan1,[scanName1 '_S0_codPPG.nii.gz']);
ppgR_S0 = niftiRead(ppgRname_S0); % correlation with PPG

ppgRname_t2s = fullfile(dDir,subj,scan1,[scanName1 '_T2s_codPPG.nii.gz']);
ppgR_t2s = niftiRead(ppgRname_t2s); % correlation with PPG


%%
%% Get ROI indices in the functional scan 
%%

% load freesurfer segmentation (these are resliced to functional space)
niFs = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_r_aseg_auto.nii.gz']));
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
% roi_list = {'CFlowvoids','AnteriorSSS','SSS','LeftTransverse','RightTransverse'};
roi_list = {'CFlowvoids','ACA','SSS','LeftTransverse','RightTransverse'};
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

    %%%% check for coordinates in functional space
    z_slice = [5 10 15 18 20 23 25 26 27 30 33 36];
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
    
    imgVol1 = ni_s0.data; %timeseries in percent signal change
    imgVol2 = ni_t2s.data;
    
    roiTS(rr).roiTrace_s0 = NaN(size(ijk_func,1),size(imgVol1,4));
    roiTS(rr).roiTrace_t2s = NaN(size(ijk_func,1),size(imgVol1,4));
    roiTS(rr).roiCod = NaN(size(ijk_func,1),2); %S0, T2s

    for kk = 1:size(ijk_func,1)
        % note that x and y are only switched around for plotting position, 
        % not for getting the actual image values
        roiTS(rr).roiTrace_s0(kk,:) = imgVol1(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3),:);
        roiTS(rr).roiTrace_t2s(kk,:) = imgVol2(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3),:);
        
        roiTS(rr).roiCod(kk,1) = ppgR_S0.data(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3));
        roiTS(rr).roiCod(kk,2) = ppgR_t2s.data(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3));

    end
end

%%%% Add freesurfer ROI traces (e.g. CSF)
fs_segm = {[3 42],[2 41],[31 63],[14],[15],[24],[4 43]};
% roi_list = {'CFlowvoids','AnteriorSSS','SSS','LeftTransverse','RightTransverse',...
%     'Gray','White','ChoroidPlexus','3rdVentr','4thVentr','CSF','LateralVentr'};
roi_list = {'CFlowvoids','ACA','SSS','LeftTransverse','RightTransverse',...
    'Gray','White','ChoroidPlexus','3rdVentr','4thVentr','CSF','LateralVentr'};

% reshape timeseries for voxel selection
tsVect1 = reshape(ni_s0.data,[numel(ni_s0.data(:,:,:,1)) length(t_me.t)]);
tsVect2 = reshape(ni_t2s.data,[numel(ni_t2s.data(:,:,:,1)) length(t_me.t)]);

% which indices are added in the output: roiTS
out_ind = 6:5+length(fs_segm);
for rr = 1:length(out_ind)
    % get the voxels from the current ROI:
    thisRoi_voxels = find(ismember(niFs.data,fs_segm{rr}));
    roiTS(out_ind(rr)).roiTrace_s0 = tsVect1(thisRoi_voxels,:);
    roiTS(out_ind(rr)).roiTrace_t2s = tsVect2(thisRoi_voxels,:);
    roiTS(out_ind(rr)).roiCod(:,1) = ppgR_S0.data(thisRoi_voxels);
    roiTS(out_ind(rr)).roiCod(:,2) = ppgR_t2s.data(thisRoi_voxels);
end
clear tsVect1 tsVect2

%% quick figure
figure('Position',[0 0 1200 400])

% remove noisy voxels from T2s here by taking roiTS(kk).roiCod(:,2)>0 ...
for kk = 1:length(roiTS)
    subplot(3,4,kk),hold on
    plot(t_me.t,mean(roiTS(kk).roiTrace_s0(roiTS(kk).roiCod(:,2)>0,:),1),'k','LineWidth',2)
    plot(t_me.t,mean(roiTS(kk).roiTrace_t2s(roiTS(kk).roiCod(:,2)>0,:),1),'Color',[.5 .5 1],'LineWidth',2)
    title(roi_list{kk})
end


%% Old code from here, we got timeseries for ROIs 
%%
%%
%%
%% save subject timeseries

tsSaveName = fullfile(dDir,subj,scan,[scanName1 '_ROItimeseries01_ME.mat']);

save(tsSaveName,'roiTS','t','roi_list')
disp(['saved ' tsSaveName])

%%
%% load group data and plot group
%%
% clear all
% close all
% roi_list = {'CFlowvoids','AnteriorSSS','SSS','LeftTransverse','RightTransverse',...
%     'Gray','White','ChoroidPlexus','3rdVentr','4thVentr','CSF','LateralVentr'};
roi_list = {'CFlowvoids','ACA','SSS','LeftTransverse','RightTransverse',...
    'Gray','White','ChoroidPlexus','3rdVentr','4thVentr','CSF','LateralVentr'};

% dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% group with FA of 48 includes sub/scan: [1/3, 2/3, 3/3, 4/1, 5/1, 6/1]

sub_nrs = [1:6];
scan_nrs = [3 3 3 1 1 1];

avg_traces = zeros(length(sub_nrs),length(roi_list),51,4);

clear out 

% load all subjects data ('roiTS','t','roi_list')
for kk = 1:length(sub_nrs)
    s_info = bb_subs(sub_nrs(kk));
    subj = s_info.subj;
    scan_nr = scan_nrs(kk);
    
    scan = s_info.scan{scan_nr};
    scanName = s_info.scanName{scan_nr};
    tsLoadName = fullfile(dDir,subj,scan,[scanName '_ROItimeseries01.mat']);

    out(kk) = load(tsLoadName,'roiTS','t','roi_list');
    
    for rr = 1:12
        avg_traces(kk,rr,:,1) = mean(out(kk).roiTS(rr).roiTrace,1);
        avg_traces(kk,rr,:,2) = mean(out(kk).roiTS(rr).roiTrace(out(kk).roiTS(rr).roiCod>.3,:),1);
        avg_traces(kk,rr,:,3) = mean(out(kk).roiTS(rr).pred,1);
        avg_traces(kk,rr,:,4) = mean(out(kk).roiTS(rr).pred(out(kk).roiTS(rr).roiCod>.3,:),1);
    end
end

%%

roi_list = {'CFlowvoids','ACA','SSS','LeftTransverse','RightTransverse',...
    'Gray','White','ChoroidPlexus','3rdVentr','4thVentr','CSF','LateralVentr'};

plot_sig = 4; % 1: mean TS, 2: mean TS for R>.3, 3: mean pred, 4: mean pred r>.3

figure('Position',[0 0 150 300])
subplot(2,1,1),hold on
plot([-1 1],[0 0],'k')
% plot Carotid flow voids (1)
roi_ind = 1;
% plot(out(1).t,squeeze(avg_traces(:,roi_ind,:,plot_sig)),'Color',[1 .5 .5],'LineWidth',1)
signal = squeeze(mean(avg_traces(:,roi_ind,:,plot_sig),1));
error = squeeze(std(avg_traces(:,roi_ind,:,plot_sig),[],1)./sqrt(length(sub_nrs)));
up_bnd = signal+error;
low_bnd = signal-error;
fill([out(1).t out(1).t(end:-1:1)],[up_bnd; low_bnd(end:-1:1)],[1 .8 .8],'EdgeColor',[1 .8 .8])
plot(out(1).t,squeeze(mean(avg_traces(:,roi_ind,:,plot_sig),1)),'r','LineWidth',2)

% plot Superior Saggital Sinus (3)
roi_ind = 3;
% plot(out(1).t,squeeze(avg_traces(:,roi_ind,:,plot_sig)),'Color',[.5 .5 1],'LineWidth',1)
signal = squeeze(mean(avg_traces(:,roi_ind,:,plot_sig),1));
error = squeeze(std(avg_traces(:,roi_ind,:,plot_sig),[],1)./sqrt(length(sub_nrs)));
up_bnd = signal+error;
low_bnd = signal-error;
fill([out(1).t out(1).t(end:-1:1)],[up_bnd; low_bnd(end:-1:1)],[.5 .5 1],'EdgeColor',[.5 .5 1])
plot(out(1).t,squeeze(mean(avg_traces(:,roi_ind,:,plot_sig),1)),'b','LineWidth',2)

xlim([-.5 .5])

subplot(2,1,2),hold on
plot([-1 1],[0 0],'k')
% % plot CSF in orange, CSF from Freesurfer has a strange location
% roi_ind = 11;
% signal = squeeze(mean(avg_traces(:,roi_ind,:,plot_sig),1));
% error = squeeze(std(avg_traces(:,roi_ind,:,plot_sig),[],1)./sqrt(length(sub_nrs)));
% up_bnd = signal+error;
% low_bnd = signal-error;
% fill([out(1).t out(1).t(end:-1:1)],[up_bnd; low_bnd(end:-1:1)],[1 .5 0],'EdgeColor',[1 1 .5])
% plot(out(1).t,squeeze(mean(avg_traces(:,roi_ind,:,plot_sig),1)),'Color',[.8 .4 0],'LineWidth',2)

% plot lateral ventricles in yellow (12)
roi_ind = 12;
% subject 2 does not have any voxels with reliable signal in this ROI...
% plot(out(1).t,squeeze(avg_traces(:,roi_ind,:,plot_sig)),'Color',[.3 .5 .3],'LineWidth',1)
signal = squeeze(mean(avg_traces(:,roi_ind,:,plot_sig),1));
error = squeeze(std(avg_traces(:,roi_ind,:,plot_sig),[],1)./sqrt(length(sub_nrs)));
up_bnd = signal+error;
low_bnd = signal-error;
fill([out(1).t out(1).t(end:-1:1)],[up_bnd; low_bnd(end:-1:1)],[1 1 .5],'EdgeColor',[1 1 .5])
plot(out(1).t,squeeze(mean(avg_traces(:,roi_ind,:,plot_sig),1)),'Color',[.8 .8 0],'LineWidth',2)


% plot choroid plexus in green (8)
roi_ind = 8; 
% plot(out(1).t,squeeze(avg_traces(:,roi_ind,:,plot_sig)),'Color',[.5 1 .5],'LineWidth',1)
signal = squeeze(mean(avg_traces(:,roi_ind,:,plot_sig),1));
error = squeeze(std(avg_traces(:,roi_ind,:,plot_sig),[],1)./sqrt(length(sub_nrs)));
up_bnd = signal+error;
low_bnd = signal-error;
fill([out(1).t out(1).t(end:-1:1)],[up_bnd; low_bnd(end:-1:1)],[.8 1 .8],'EdgeColor',[.8 1 .8])
plot(out(1).t,squeeze(mean(avg_traces(:,roi_ind,:,plot_sig),1)),'g','LineWidth',2)
xlim([-.5 .5])

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','ROI',['TimeSerie_6subAvg-FA48_meanTS_COD0_3']))
% print('-painters','-r300','-depsc',fullfile(dDir,'figures','ROI',['TimeSerie_6subAvg-FA48_meanTS_COD0_3']))
% print('-painters','-r300','-dpng',[dDir './figures/ROI/TimeSerie_6subAvg-FA48_meanTS'])
% print('-painters','-r300','-depsc',[dDir './figures/ROI/TimeSerie_6subAvg-FA48_meanTS'])

print('-painters','-r300','-dpng',fullfile(dDir,'figures','ROI',['TimeSerie_6subAvg-FA48_meanPRED_COD0_3']))
print('-painters','-r300','-depsc',fullfile(dDir,'figures','ROI',['TimeSerie_6subAvg-FA48_meanPRED_COD0_3']))

%%
figure
for kk = 1:length(sub_nrs)
    for rr = 1:12
        subplot(6,2,rr),hold on
        plot(out(kk).t,mean(out(kk).roiTS(rr).roiTrace,1))
        % only voxels with R>.3
%         plot(out(kk).t,mean(out(kk).roiTS(rr).roiTrace(out(kk).roiTS(rr).roiCod>.3,:),1))
    end
end

for rr = 1:12
    subplot(6,2,rr)
    ylabel(out(1).roi_list{rr})
end

%% plot pc1 weight
figure
for kk = 1:length(sub_nrs)
    for rr = 1:11
        subplot(6,2,rr),hold on
        plot(kk,median(out(kk).roiTS(rr).roiWeights(out(kk).roiTS(rr).roiCod>.3,1),1),'.')
%         plot(kk,out(kk).roiTS(rr).roiWeights(:,1),'.')
%        boxplot(out(kk).roiTS(rr).roiWeights(:,1))
    end
end

for rr = 1:11
    subplot(6,2,rr)
    plot([0 7],[0 0],'k')
    ylabel(out(1).roi_list{rr})
end


% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/ROI/TimeSerie_sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_roi-CSF_COD-0_3'])
% print('-painters','-r300','-depsc',[dDir './figures/ROI/TimeSerie_sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_roi-CSF_COD-0_3'])

