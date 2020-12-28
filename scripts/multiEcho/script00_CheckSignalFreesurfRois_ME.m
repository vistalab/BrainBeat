clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';
% chdir(dDir)

%% Generate basic fMRI change from Freesurfer segmentations 

% Functionals
% Select a subject and scan nummer
s_nr = 4; %[4 5 6]

% get multi echo data
if s_nr==4
    scan_nrs = {[4 5],[6 7]};
elseif s_nr==5
    scan_nrs = {[4 5],[6 7]};
elseif s_nr==6
    scan_nrs = {[4 5],[6 7],[9 10]};
end

s_info = bb_subs(s_nr);
subj = s_info.subj;

% Get the anatomical:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));


scn_me = 1;%:length(scan_nrs)
% first echo
scan1 = s_info.scan{scan_nrs{scn_me}(1)};
scanName1 = s_info.scanName{scan_nrs{scn_me}(1)};

% load first dataset to get PPG curve here:
fmri = fullfile(dDir,subj,scan1,[scanName1 '.nii.gz']);
if ~exist(fmri,'file')
    clear ni
    error('filename %s does not exist',fmri)
end
ni = niftiRead(fmri);

% load coregistration matrix:
load(fullfile(dDir,subj,scan1,[scanName1 'AcpcXform_new.mat']));
acpcXform = acpcXform_new;

% load time T ppgT/t_me:
t_me = load(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponseT']),'t');

% load average of all heartbeats:
fname_s0 = [scanName1 '_PPGtrigResponse_S0.nii.gz'];
fname_t2s = [scanName1 '_PPGtrigResponse_T2s.nii.gz'];
ni_s0 = niftiRead(fullfile(dDir,subj,scan1,fname_s0));
ni_t2s = niftiRead(fullfile(dDir,subj,scan1,fname_t2s));

% load the correlation with heartbeat (made with bbCorrelate2physio):
%%% DO WE NEED THIS?
% ppgRname_S0 = fullfile(dDir,subj,scan1,[scanName1 '_S0_codPPG.nii.gz']);
% ppgR_S0 = niftiRead(ppgRname_S0); % correlation with PPG
% ppgRname_t2s = fullfile(dDir,subj,scan1,[scanName1 '_T2s_codPPG.nii.gz']);
% ppgR_t2s = niftiRead(ppgRname_t2s); % correlation with PPG

% get physio stuff we need:
physio      = physioCreate('nifti',ni);
ppgCurve    = physioGet(physio,'ppg ppgcurve');
ppgCurveT   = physioGet(physio,'ppg ppgtcurve');
srate       = 1/bbGet(ni,'tr');

niSegm = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_r_DKTatlas_aseg.nii.gz']));
% matrix to vector
segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

roiNames = {'lingual','insula','central','cingulate'};
roiCodes = {[1013 2013],[1035 2035],[1022 1024 2022 2024],[1026 2026]};

% put the data in a matrix Voxel X Time
respMat_s0 = reshape(ni_s0.data,[size(ni_s0.data,1) * size(ni_s0.data,2) * size(ni_s0.data,3)],size(ni_s0.data,4));
respMat_t2s = reshape(ni_t2s.data,[size(ni_t2s.data,1) * size(ni_t2s.data,2) * size(ni_t2s.data,3)],size(ni_t2s.data,4));

avResp_s0 = zeros(size(ni_s0.data,4),length(roiNames));
avResp_t2s = zeros(size(ni_t2s.data,4),length(roiNames));
for kk = 1:length(roiNames)
    avResp_s0(:,kk) = mean(respMat_s0(ismember(segmVect,roiCodes{kk}),:),1);
    avResp_t2s(:,kk) = mean(respMat_t2s(ismember(segmVect,roiCodes{kk}),:),1);
end


%%
areas_plot = [1:4];

figure('Position',[0 0 300 600])

for kk = 1:length(areas_plot)
    thisArea = areas_plot(kk);
    subplot(length(areas_plot)+1,1,kk),hold on
    plot(t_me.t,avResp_s0(:,thisArea)*100,'Color',[0 .6 .8],'LineWidth',2)
    plot(t_me.t,avResp_t2s(:,thisArea)*100,'Color',[.3 .8 .6],'LineWidth',2)
%     plot([0 0],[min(avResp(:,thisArea)*100) max(avResp(:,thisArea)*100)],'k')
    xlim([min(t_me.t) max(t_me.t)])
    title(roiNames{thisArea})
end

% ppg for completeness
subplot(length(areas_plot)+1,1,length(areas_plot)+1),hold on
plot(ppgCurveT,ppgCurve,'Color',[0 .6 .8],'LineWidth',2)
plot([0 0],[min(ppgCurve) max(ppgCurve)],'k')
xlim([min(ppgCurveT) max(ppgCurveT)])
title('ppg')

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir '/figures/segmentation/s' int2str(s_nr) '_scan' int2str(scan_nr) '_tracessegm'])
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/s' int2str(s_nr) '_scan' int2str(scan_nr) '_tracessegm'])




%% Now plot the average responses across subjects

all_subs = [4 5 6];

% roiNames = {'lingual','insula','anteriorcingulate','central','middlefrontal','itg','fusiform'};
% roiCodes = {[1013 2013],[1035 2035],[1026 2026],[1022 1024 2022 2024],[1027 2027],[1011 2011],[1007 2007]};
% roiNames = {'3rd Ventr','4th Ventr','lateral ventr','inferior lat ventr'};
% roiCodes = {[14],[15],[4 43],[5 44]};
% get responses for these ROIs
roiNames = {'lingual','insula','cingulate','central'};
roiCodes = {[1013 2013],[1035 2035],[1026 2026],[1022 1024 2022 2024]};

t_hr = linspace(-.5,1.5,128);
avResp_hr_s0 = NaN(length(t_hr),length(roiNames),length(all_subs));
avResp_hr_t2s = NaN(length(t_hr),length(roiNames),length(all_subs));
ppgCurve_hr = NaN(length(t_hr),length(all_subs));

    
for ss = 1:length(all_subs)
    % Select a subject and scan nummer
    s_nr = all_subs(ss); %[4 5 6]
    
    % get multi echo data
    if s_nr==4
        scan_nrs = {[4 5],[6 7]};
    elseif s_nr==5
        scan_nrs = {[4 5],[6 7]};
    elseif s_nr==6
        scan_nrs = {[4 5],[6 7],[9 10]};
    end

    s_info = bb_subs(s_nr);
    subj = s_info.subj;

    % Get the anatomical:
    niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

    scn_me = 1;%:length(scan_nrs)
    % first echo
    scan1 = s_info.scan{scan_nrs{scn_me}(1)};
    scanName1 = s_info.scanName{scan_nrs{scn_me}(1)};

    % load first dataset to get PPG curve here:
    fmri = fullfile(dDir,subj,scan1,[scanName1 '.nii.gz']);
    if ~exist(fmri,'file')
        clear ni
        error('filename %s does not exist',fmri)
    end
    ni = niftiRead(fmri);

    % get physio stuff we need:
    physio      = physioCreate('nifti',ni);
    ppgCurve    = physioGet(physio,'ppg ppgcurve');
    ppgCurveT   = physioGet(physio,'ppg ppgtcurve');
    srate       = 1/bbGet(ni,'tr');

    % load time T ppgT/t_me:
    t_me = load(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponseT']),'t');

    % load average of all heartbeats:
    fname_s0 = [scanName1 '_PPGtrigResponse_S0.nii.gz'];
    fname_t2s = [scanName1 '_PPGtrigResponse_T2s.nii.gz'];
    ni_s0 = niftiRead(fullfile(dDir,subj,scan1,fname_s0));
    ni_t2s = niftiRead(fullfile(dDir,subj,scan1,fname_t2s));
    
    % load segmentation
    niSegm = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_r_DKTatlas_aseg.nii.gz']));
    % matrix to vector
    segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

    % put the data in a matrix Voxel X Time
    respMat_s0 = reshape(ni_s0.data,[size(ni_s0.data,1) * size(ni_s0.data,2) * size(ni_s0.data,3)],size(ni_s0.data,4));
    respMat_t2s = reshape(ni_t2s.data,[size(ni_t2s.data,1) * size(ni_t2s.data,2) * size(ni_t2s.data,3)],size(ni_t2s.data,4));

    avResp_s0 = zeros(size(ni_s0.data,4),length(roiNames));
    avResp_t2s = zeros(size(ni_t2s.data,4),length(roiNames));
    for kk = 1:length(roiNames)
        avResp_s0(:,kk) = mean(respMat_s0(ismember(segmVect,roiCodes{kk}),:),1);
        avResp_t2s(:,kk) = mean(respMat_t2s(ismember(segmVect,roiCodes{kk}),:),1);
    end

    % resample avResp to ppg cycle
    ppg_cycle   = 1./physioGet(physio,'PPGrate');
    avResp_sel_s0 = avResp_s0(t_me.t>=(0-(.5*ppg_cycle)) & t_me.t<=1.5*ppg_cycle,:);
    avResp_sel_t2s = avResp_t2s(t_me.t>=(0-(.5*ppg_cycle)) & t_me.t<=1.5*ppg_cycle,:);
    t_sel = t_me.t(t_me.t>=(0-(.5*ppg_cycle)) & t_me.t<=1.5*ppg_cycle);
    ppgCurve_sel = ppgCurve(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);
    ppgt_sel = ppgCurveT(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);

    t_temp = linspace(min(t_sel),max(t_sel),128);

    avResp_hr_s0(:,:,ss) = interp1(t_sel,avResp_sel_s0,t_temp);
    avResp_hr_t2s(:,:,ss) = interp1(t_sel,avResp_sel_t2s,t_temp);
    ppgCurve_hr(:,ss) = interp1(ppgt_sel,ppgCurve_sel,t_temp);

    clear ni_s0 ni_t2s t_me
end


%%
areas_plot = [1:4];
figure('Position',[0 0 240 700])
for kk = 1:length(areas_plot)
    thisArea = areas_plot(kk);
    subplot(length(areas_plot)+1,2,kk*2-1),hold on
    plot(t_hr,squeeze(avResp_hr_s0(:,thisArea,:))*100,'Color',[0 .6 .8])
    plot(t_hr,median(avResp_hr_s0(:,thisArea,:),3)*100,'Color',[0 .6 .8],'LineWidth',2)

    plot([0 0],[min(mean(avResp_hr_s0(:,thisArea,:),3)*100)*2 max(mean(avResp_hr_s0(:,thisArea,:),3)*100)*2],'k')
    xlim([min(t_hr) max(t_hr)])
    ylabel(roiNames{thisArea})
    
    subplot(length(areas_plot)+1,2,kk*2),hold on
    plot(t_hr,squeeze(avResp_hr_t2s(:,thisArea,:))*100,'Color',[.3 .8 .6])
    plot(t_hr,median(avResp_hr_t2s(:,thisArea,:),3)*100,'Color',[.3 .8 .6],'LineWidth',2)

    plot([0 0],[min(mean(avResp_hr_s0(:,thisArea,:),3)*100)*2 max(mean(avResp_hr_s0(:,thisArea,:),3)*100)*2],'k')
    xlim([min(t_hr) max(t_hr)])
    
%     ylim([-.5 .5])
end

subplot(length(areas_plot)+1,2,2*length(areas_plot)+1),hold on
plot(t_hr,ppgCurve_hr./nanstd(ppgCurve_hr,[],1),'Color',[0 .6 .8])
plot(t_hr,mean(ppgCurve_hr./nanstd(ppgCurve_hr,[],1),2),'Color',[0 .6 .8],'LineWidth',2)
plot([0 0],[-2 4],'k')
xlim([min(t_hr) max(t_hr)])
ylabel('PPG')
subplot(length(areas_plot)+1,2,2*length(areas_plot)+2),hold on
plot(t_hr,ppgCurve_hr./nanstd(ppgCurve_hr,[],1),'Color',[0 .6 .8])
plot(t_hr,mean(ppgCurve_hr./nanstd(ppgCurve_hr,[],1),2),'Color',[0 .6 .8],'LineWidth',2)
plot([0 0],[-2 4],'k')
xlim([min(t_hr) max(t_hr)])


% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir '/figures/segmentation/allsubjects_FA48_tracesDKT1'])
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/allsubjects_FA48_tracesDKT1'])

%%
%% Now plot the average responses across subjects

all_subs = [4 5 6];

t_hr = linspace(-.5,1.5,128);
avResp_hr = NaN(length(t_hr),5,length(all_subs));
ppgCurve_hr = NaN(length(t_hr),length(all_subs));
roiNames = {'GM','WM','Ventricles','CSF','Veno'};
roiCodes = {[1],[2],[3],[4],[5]};

t_hr = linspace(-.5,1.5,128);
avResp_hr_s0 = NaN(length(t_hr),length(roiNames),length(all_subs));
avResp_hr_t2s = NaN(length(t_hr),length(roiNames),length(all_subs));

for ss = 1:length(all_subs)
    % Select a subject and scan nummer
    s_nr = all_subs(ss); %[4 5 6]
    
    % get multi echo data
    if s_nr==4
        scan_nrs = {[4 5],[6 7]};
    elseif s_nr==5
        scan_nrs = {[4 5],[6 7]};
    elseif s_nr==6
        scan_nrs = {[4 5],[6 7],[9 10]};
    end
    
    s_info = bb_subs(s_nr);
    subj = s_info.subj;

    % Get the anatomical:
    niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

    scn_me = 1;%:length(scan_nrs)
    % first echo
    scan1 = s_info.scan{scan_nrs{scn_me}(1)};
    scanName1 = s_info.scanName{scan_nrs{scn_me}(1)};

    % load first dataset to get PPG curve here:
    fmri = fullfile(dDir,subj,scan1,[scanName1 '.nii.gz']);
    if ~exist(fmri,'file')
        clear ni
        error('filename %s does not exist',fmri)
    end
    ni = niftiRead(fmri);

    % get physio stuff we need:
    physio      = physioCreate('nifti',ni);
    ppgCurve    = physioGet(physio,'ppg ppgcurve');
    ppgCurveT   = physioGet(physio,'ppg ppgtcurve');
    srate       = 1/bbGet(ni,'tr');

    % load time T ppgT/t_me:
    t_me = load(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponseT']),'t');

    % load average of all heartbeats:
    fname_s0 = [scanName1 '_PPGtrigResponse_S0.nii.gz'];
    fname_t2s = [scanName1 '_PPGtrigResponse_T2s.nii.gz'];
    ni_s0 = niftiRead(fullfile(dDir,subj,scan1,fname_s0));
    ni_t2s = niftiRead(fullfile(dDir,subj,scan1,fname_t2s));
    
    % load segmentation: Freesurfer for GM, WM, Ventricles, CSF from SPM and Venogram
    niSegm = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_combineSegm.nii.gz']));   
    % matrix to vector
    segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

    % put the data in a matrix Voxel X Time
    respMat_s0 = reshape(ni_s0.data,[size(ni_s0.data,1) * size(ni_s0.data,2) * size(ni_s0.data,3)],size(ni_s0.data,4));
    respMat_t2s = reshape(ni_t2s.data,[size(ni_t2s.data,1) * size(ni_t2s.data,2) * size(ni_t2s.data,3)],size(ni_t2s.data,4));

    avResp_s0 = zeros(size(ni_s0.data,4),length(roiNames));
    avResp_t2s = zeros(size(ni_t2s.data,4),length(roiNames));
    for kk = 1:length(roiNames)
%         avResp_s0(:,kk) = mean(respMat_s0(ismember(segmVect,roiCodes{kk}),:),1);
%         avResp_t2s(:,kk) = mean(respMat_t2s(ismember(segmVect,roiCodes{kk}),:),1);
        
        % try with resample
        temp_s0 = NaN(size(ni_s0.data,4),100);
        temp_t2s = NaN(size(ni_t2s.data,4),100);
        for n_rep = 1:1000
            a_s0 = respMat_s0(ismember(segmVect,roiCodes{kk}),:);
            a_t2s = respMat_t2s(ismember(segmVect,roiCodes{kk}),:);
            b = randsample(1:size(a_s0,1),round(size(a_s0,1)/4),true);
            temp_s0(:,n_rep) = mean(a_s0(b,:),1);
            temp_t2s(:,n_rep) = mean(a_t2s(b,:),1);
        end
        avResp_s0(:,kk) = mean(temp_s0,2);
        avResp_t2s(:,kk) = mean(temp_t2s,2);
    end

    % resample avResp to ppg cycle
    ppg_cycle   = 1./physioGet(physio,'PPGrate');
    avResp_sel_s0 = avResp_s0(t_me.t>=(0-(.5*ppg_cycle)) & t_me.t<=1.5*ppg_cycle,:);
    avResp_sel_t2s = avResp_t2s(t_me.t>=(0-(.5*ppg_cycle)) & t_me.t<=1.5*ppg_cycle,:);
    t_sel = t_me.t(t_me.t>=(0-(.5*ppg_cycle)) & t_me.t<=1.5*ppg_cycle);
    ppgCurve_sel = ppgCurve(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);
    ppgt_sel = ppgCurveT(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);

    t_temp = linspace(min(t_sel),max(t_sel),128);

    avResp_hr_s0(:,:,ss) = interp1(t_sel,avResp_sel_s0,t_temp);
    avResp_hr_t2s(:,:,ss) = interp1(t_sel,avResp_sel_t2s,t_temp);
    ppgCurve_hr(:,ss) = interp1(ppgt_sel,ppgCurve_sel,t_temp);

    clear ni_s0 ni_t2s t_me
end

%%
areas_plot = [1 3 5];
figure('Position',[0 0 240 400])
for kk = 1:length(areas_plot)
    thisArea = areas_plot(kk);
    subplot(length(areas_plot)+1,2,kk*2-1),hold on
    plot(t_hr,squeeze(avResp_hr_s0(:,thisArea,:))*100,'Color',[0 .6 .8])
    plot(t_hr,mean(avResp_hr_s0(:,thisArea,:),3)*100,'Color',[0 .6 .8],'LineWidth',2)
    plot([0 0],[min(mean(avResp_hr_s0(:,thisArea,:),3)*100)*2 max(mean(avResp_hr_s0(:,thisArea,:),3)*100)*2],'k')
    xlim([min(t_hr) max(t_hr)])
    ylabel(roiNames{thisArea})

    subplot(length(areas_plot)+1,2,kk*2),hold on
    plot(t_hr,squeeze(avResp_hr_t2s(:,thisArea,:))*100,'Color',[0 .6 .8])
    plot(t_hr,mean(avResp_hr_t2s(:,thisArea,:),3)*100,'Color',[0 .6 .8],'LineWidth',2)
    plot([0 0],[min(mean(avResp_hr_t2s(:,thisArea,:),3)*100)*2 max(mean(avResp_hr_t2s(:,thisArea,:),3)*100)*2],'k')
    xlim([min(t_hr) max(t_hr)])
end

subplot(length(areas_plot)+1,2,2*length(areas_plot)+1),hold on
plot(t_hr,ppgCurve_hr./nanstd(ppgCurve_hr,[],1),'Color',[0 .6 .8])
plot(t_hr,mean(ppgCurve_hr./nanstd(ppgCurve_hr,[],1),2),'Color',[0 .6 .8],'LineWidth',2)
plot([0 0],[-2 4],'k')
xlim([min(t_hr) max(t_hr)])
ylabel('PPG')

subplot(length(areas_plot)+1,2,2*length(areas_plot)+2),hold on
plot(t_hr,ppgCurve_hr./nanstd(ppgCurve_hr,[],1),'Color',[0 .6 .8])
plot(t_hr,mean(ppgCurve_hr./nanstd(ppgCurve_hr,[],1),2),'Color',[0 .6 .8],'LineWidth',2)
plot([0 0],[-2 4],'k')
xlim([min(t_hr) max(t_hr)])


% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir '/figures/segmentation/allsubjects_FA48_tracessegm'])
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/allsubjects_FA48_tracessegm'])

