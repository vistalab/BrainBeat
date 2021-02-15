clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';
% chdir(dDir)

%% Basic check from here!!!
%% Generate basic fMRI change from segmentation 

% Functionals
% Select a subject and scan nummer
s_nr = 1; %[1 2 3 4 5 6]
scan_nr = 3; % [3 3 3 1 1 1]

dkt_table = readtable('dkt_areas.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dktNames = dkt_table.label;
dktCodes = dkt_table.label_nr;

subs = bb_subs(s_nr);
subj = subs.subj;
scan = subs.scan{scan_nr};
scanName = subs.scanName{scan_nr};
fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
if ~exist(fmri,'file')
    clear ni
    error('filename %s does not exist',fmri)
end
ni = niftiRead(fmri);

% get physio stuff we need:
physio      = physioCreate('nifti',ni);
ppg_onsets  = physioGet(physio,'ppg peaks');
ppgCurve    = physioGet(physio,'ppg ppgcurve');
ppgCurveT   = physioGet(physio,'ppg ppgtcurve');
srate       = 1/bbGet(ni,'tr');

ppgResp = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse.nii.gz']));
ppgT = load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT.mat']));

niSegm = niftiRead(fullfile(dDir,subj,scan,[scanName '_combineSegm.nii.gz']));
roiNames = {'GM','WM','Ventricles','CSF','Veno'};
% Freesurfer for GM, WM, Ventricles, CSF from SPM and Venogram

% Voxel vector with segmentation labels 1:5 for [GM WM Ventricles CSF Veno]
segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

% load DKT atlas
niDKT = niftiRead(fullfile(dDir,subj,scan,[scanName '_r_DKTatlas_aseg.nii.gz']));
% matrix to vector
segmVectDKT = reshape(niDKT.data,[size(niDKT.data,1) * size(niDKT.data,2) * size(niDKT.data,3)],1);

%%%% put the PPG triggered data in a matrix Voxel X Time
respMat = reshape(ppgResp.data,[size(ppgResp.data,1) * size(ppgResp.data,2) * size(ppgResp.data,3)],size(ppgResp.data,4));

% initialize average PPG triggered response across ROIs
avResp = zeros(size(ppgResp.data,4),length(roiNames)+1);
for kk = 1:length(roiNames)
    avResp(:,kk) = mean(respMat(segmVect==kk,:),1);
end

this_areas = 1; % cingulate
avResp(:,kk+1) = mean(respMat(ismember(segmVectDKT,dktCodes(this_areas)),:),1);


%%%% put the fMRI data in a matrix Voxel X Time
sigMat = reshape(ni.data,[size(ni.data,1) * size(ni.data,2) * size(ni.data,3)],size(ni.data,4));

% initialize average total time signal
avSig = zeros(size(ni.data,4),length(roiNames)+1);
for kk = 1:length(roiNames)
    avSig(:,kk) = mean(sigMat(segmVect==kk,:),1);
end
avSig(:,kk+1) = mean(sigMat(ismember(segmVectDKT,dktCodes(this_areas)),:),1);
avSig(1:5,:) = NaN; % first scans to NaN

roiNames(length(roiNames)+1) = dktNames(this_areas);

%%

areas_plot = [6 3 5];

tt = [1:size(avSig,1)]/srate;

figure('Position',[0 0 500 400])
for kk = 1:length(areas_plot)
    thisArea = areas_plot(kk);
    
    sig_plot = avSig(:,thisArea)';
    % detrend
    x = 6:length(sig_plot);
    p = polyfit(x,sig_plot(x),1);    

    % percent modulation
    mean_factor = mean(sig_plot(x)); % mean
    sig_plot = 100*(sig_plot - (p(1)*[1:length(sig_plot)] + p(2)))./mean_factor;
    
    subplot(length(areas_plot)+1,2,kk*2-1),hold on
    plot(tt,sig_plot,'.','Color',[0 .6 .8],'MarkerSize',10)
    plot(tt,sig_plot,'Color',[0 .6 .8],'LineWidth',1)
%     plot([ppg_onsets ppg_onsets],[nanmean(avSig(:,thisArea))*100 nanmean(avSig(:,thisArea))*100],'k')
    plot([ppg_onsets ppg_onsets],[min(avResp(:,thisArea)*100) max(avResp(:,thisArea)*100)],'k')
    xlim([85 95]) 
end

for kk = 1:length(areas_plot)
    thisArea = areas_plot(kk);
    subplot(length(areas_plot)+1,2,kk*2),hold on
    plot(ppgT.t,avResp(:,thisArea)*100,'Color',[0 .6 .8],'LineWidth',2)
    plot([0 0],[min(avResp(:,thisArea)*100) max(avResp(:,thisArea)*100)],'k')
    xlim([min(ppgT.t) max(ppgT.t)])
    title(roiNames{thisArea})
end

% ppg for completeness
subplot(length(areas_plot)+1,2,length(areas_plot)*2+1),hold on
plot([ppg_onsets ppg_onsets],[-50 100],'k')
plot([1:length(physio.ppg.data)]/physio.ppg.srate,physio.ppg.data)
xlim([85 95])

subplot(length(areas_plot)+1,2,length(areas_plot)*2+2),hold on
plot(ppgCurveT,ppgCurve,'Color',[0 .6 .8],'LineWidth',2)
plot([0 0],[min(ppgCurve) max(ppgCurve)],'k')
xlim([min(ppgCurveT) max(ppgCurveT)])
title('ppg')

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir '/figures/segmentation/s' int2str(s_nr) '_scan' int2str(scan_nr) '_tracessegm'])
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/s' int2str(s_nr) '_scan' int2str(scan_nr) '_tracessegm'])


%%

figure,
for kk = 1:5
    subplot(1,5,kk)
    imagesc(ppgT.t,[1:length(find(segmVect==kk))],100*respMat(segmVect==kk,:),[-1 1])
    hold on
    plot([0 0],[1 length(find(segmVect==kk))],'k','LineWidth',2)
end


%% Now plot the average responses across subjects

all_subs = [1 2 3 4 5 6];
all_scans = [3 3 3 1 1 1];

t_hr = linspace(-.5,1.5,128);
avResp_hr = NaN(length(t_hr),5,length(all_subs));
ppgCurve_hr = NaN(length(t_hr),length(all_subs));
roiNames = {'GM','WM','Ventricles','CSF','Veno'};
    
for ss = 1:length(all_subs)
    % Functionals
    % Select a subject and scan nummer
    s_nr = all_subs(ss); %[1 2 3 4 5 6]
    scan_nr = all_scans(ss); % [3 3 3 1 1 1]

    subs = bb_subs(s_nr);
    subj = subs.subj;
    scan = subs.scan{scan_nr};
    scanName = subs.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
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

    ppgResp = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse.nii.gz']));
    ppgT = load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT.mat']));

    niSegm = niftiRead(fullfile(dDir,subj,scan,[scanName '_combineSegm.nii.gz']));
    % Freesurfer for GM, WM, Ventricels, CSF from SPM and Venogram

    % Voxel vector with segmentation labels 1:5 for [SM WM Ventricles CSF Veno]
    segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

    % put the data in a matrix Voxel X Time
    respMat = reshape(ppgResp.data,[size(ppgResp.data,1) * size(ppgResp.data,2) * size(ppgResp.data,3)],size(ppgResp.data,4));

    avResp = zeros(size(ppgResp.data,4),length(roiNames));
    for kk = 1:length(roiNames)
        avResp(:,kk) = mean(respMat(segmVect==kk,:),1);
    end

    % put the data in a matrix Voxel X Time
    sigMat = reshape(ni.data,[size(ni.data,1) * size(ni.data,2) * size(ni.data,3)],size(ni.data,4));

    avSig = zeros(size(ni.data,4),length(roiNames));
    for kk = 1:length(roiNames)
        avSig(:,kk) = mean(sigMat(segmVect==kk,:),1);
    end
    avSig(1:5,:) = NaN; % first scans to NaN

    % resample avResp to ppg cycle
    ppg_cycle   = 1./physioGet(physio,'PPGrate');
    avResp_sel = avResp(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    t_sel = ppgT.t(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle);
    ppgCurve_sel = ppgCurve(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);
    ppgt_sel = ppgCurveT(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);

    t_temp = linspace(min(t_sel),max(t_sel),128);

    avResp_hr(:,:,ss) = interp1(t_sel,avResp_sel,t_temp);
    ppgCurve_hr(:,ss) = interp1(ppgt_sel,ppgCurve_sel,t_temp);

end


%%
areas_plot = [1 3 5];
figure('Position',[0 0 120 400])
for kk = 1:length(areas_plot)
    thisArea = areas_plot(kk);
    subplot(length(areas_plot)+1,1,kk),hold on
    plot(t_hr,squeeze(avResp_hr(:,thisArea,:))*100,'Color',[0 .6 .8])
    plot(t_hr,mean(avResp_hr(:,thisArea,:),3)*100,'Color',[0 .6 .8],'LineWidth',2)
    plot([0 0],[min(mean(avResp_hr(:,thisArea,:),3)*100)*2 max(mean(avResp_hr(:,thisArea,:),3)*100)*2],'k')
    xlim([min(t_hr) max(t_hr)])
    title(roiNames{thisArea})
end

subplot(length(areas_plot)+1,1,length(areas_plot)+1),hold on
plot(t_hr,ppgCurve_hr./nanstd(ppgCurve_hr,[],1),'Color',[0 .6 .8])
plot(t_hr,mean(ppgCurve_hr./nanstd(ppgCurve_hr,[],1),2),'Color',[0 .6 .8],'LineWidth',2)
plot([0 0],[-2 4],'k')
xlim([min(t_hr) max(t_hr)])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir '/figures/segmentation/allsubjects_FA48_tracessegm'])
print('-painters','-r300','-depsc',[dDir '/figures/segmentation/allsubjects_FA48_tracessegm'])


%%
%%
%%
%% Other checks:
%%
%% The T2* data are here.  

% Select a subject and scan nummer
s_nr = 1;
scan_nr = 3;

subs = bb_subs(s_nr);
subj = subs.subj;
scan = subs.scan{scan_nr};
scanName = subs.scanName{scan_nr};
fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
if ~exist(fmri,'file')
    clear ni
    error('filename %s does not exist',fmri)
end
ni = niftiRead(fmri);

% load coregistration matrix for the functionals:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

% Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,subs.veno,[subs.venoName '.nii']));
% % load coregistration matrix (for the venogram):
% xf_veno=load(fullfile(dDir,subj,subs.veno,[subs.venoName 'AcpcXform.mat']));

%% Anatomicals

anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii.gz']);
niAnatomy = niftiRead(anat);

%% Overlay functionals and anatomy to check coregistration

sliceThisDim = 3;
if s_nr == 1
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [1,10,-20];
%     curPos = [-29,-14,-50]; % 03
%     curPos = [-3,30,-43]; % 03
    curPos = [1,10,-20];
    curPos = [1,1,-20];
%     curPos = [-11 34 -71]; % Carotid
elseif s_nr == 2
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
%     curPos = [0,4,38];
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
end
niFunc = ni;
physio = physioCreate('nifti',ni);

% niFunc.data = ni.data(:,:,:,1); % overlay the first functional - more structure visible
% bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
% title('First functional on anatomy')
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/checkCoreg/' subj '_' scan '_Func1onAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

niFunc.data = mean(ni.data(:,:,:,5:end),4); % overlay the mean functional
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
title('Mean functional on anatomy')
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/checkCoreg/' subj '_' scan '_MeanFuncOnAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

ppgTSname = fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT.mat']),'t');

%%
all_xyz = {[27 51 20],[27 50 20],[28 51 20],[28 50 20]};

figure
subplot(2,1,1)
imagesc(niFunc.data(:,:,20)),hold on

for kk = 1:length(all_xyz)
    xyz = all_xyz{kk};
    subplot(2,1,1)
    plot(xyz(2),xyz(1),'k.')
    axis image

    subplot(2,1,2),hold on
    srate = 1/bbGet(ni,'tr');
    ttt = (1:size(ni.data,4))/srate;
    thisSignal = squeeze(ni.data(xyz(1),xyz(2),xyz(3),:));
    plot(ttt,thisSignal)
end


ppg_onsets = physioGet(physio,'ppg peaks');
plot(ppg_onsets,mean(thisSignal),'r.')

% get scans nrs closest to PPG peak, get diff
ppg_scannr_onsets = round(ppg_onsets*srate);
sig_slope = diff(ni.data,[],4);
%%

figure
subplot(1,3,1)
slice_nr = 17;
slice_mean_df = mean(sig_slope(:,:,slice_nr,ppg_scannr_onsets),4);
imagesc(slice_mean_df,[-3 3])
subplot(1,3,2)
slice_nr = 29;
slice_mean_df = squeeze(mean(sig_slope(slice_nr,:,:,ppg_scannr_onsets),4));
imagesc(slice_mean_df',[-3 3]),axis xy

%% on anatomical


sliceThisDim = 3;
imDims = [-90 -120 -120; 90 130 90];
curPos = [1,1,-20];
niDF = niFunc;
niDF.data = mean(sig_slope,4);
niDF.data = -min(niDF.data(:))+niDF.data;
bbOverlayFuncAnat(niDF,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
title('Mean change direction at PPG peak on anatomy')
set(gcf,'PaperPositionMode','auto')


sliceThisDim = 1;
imDims = [-90 -120 -120; 90 130 90];
curPos = [1,1,-20];
niDF = niFunc;
niDF.data = mean(sig_slope,4);
niDF.data = -min(niDF.data(:))+niDF.data;
bbOverlayFuncAnat(niDF,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.35)
title('Mean change direction at PPG peak on anatomy')
set(gcf,'PaperPositionMode','auto')

%% 

exampl_coords = [27 51 20];%[27 51 20] = superior saggital sinus, S1, scan3
% exampl_coords = [27 45 26];%[27 51 20] = superior saggital sinus, S2, scan3
plot_nrs = [13:2:20]; % heartbeat nrs to plot

figure('Position',[0 0 250 400])

slices=[1:size(ni.data,3)];

% get the nifti stuff we need:
timing=bbGet(ni,'timing');
mux_f=bbGet(ni,'super slices');
srate=1/bbGet(ni,'tr');

% get physio stuff we need:
physio     = physioCreate('nifti',ni);
ppg_onsets = physioGet(physio,'ppg peaks');

% epoch definitions
epoch_pre = .5;%sec pre-onset
epoch_post = 1;%sec post-onset

step_size = 1/srate/mux_f;% in s
srate_epochs = 1/step_size;
t = [-epoch_pre:step_size:epoch_post];%s timing for 1 epoch
t_vox = min(timing(:)):step_size:max(timing(:)); % maximal timing accuracy for this scan session

% only use those ppg onsets that have an entire trial before and after
ppg_onsets = ppg_onsets((ppg_onsets-epoch_pre)>1); % get rid of early ones, and the first scans
ppg_onsets = ppg_onsets((ppg_onsets+epoch_post)<max(timing(:))); % get rid of late ones

for s = exampl_coords(3)%1:length(slices)
    disp(['slice ' int2str(s) ' of ' int2str(length(slices))])
    sli = slices(s);
    d = squeeze(ni.data(:,:,sli,1:end));

    % demean and express in percent modulation 
    d_norm=reshape(d,[size(d,1) * size(d,2), size(d,3)]);
    points_use=4:size(d_norm,2); % do not use the first couple of scans
    for kk = 1:size(d_norm,1)
        % do not use first scans:
        x = points_use;
        y = d_norm(kk,points_use);
        
        % detrend
        p = polyfit(x,y,1);    

        % percent modulation
        mean_factor = mean(d_norm(kk,points_use)); % mean
        d_norm(kk,:) = (d_norm(kk,:) - (p(1)*[1:size(d_norm,2)] + p(2)))./mean_factor;
    end
    d = 100*reshape(d_norm,[size(d,1), size(d,2), size(d,3)]);
    clear d_norm

    % get the timing for this slice
    t_sli = timing(sli,:);

    subplot(4,1,1),hold on
    plot([1:length(physio.ppg.data)]./physio.ppg.srate,physio.ppg.data,'k')
    %  plot(ppg_onsets,physio.ppg.data(round(ppg_onsets*physio.ppg.srate)),'ko','MarkerSize',10)
    plot([ppg_onsets ppg_onsets],[-20 50],'k:') 
    xlim([ppg_onsets(plot_nrs(1))-epoch_pre ppg_onsets(plot_nrs(end)+1)+epoch_post])
    title('PPG signal')
%     ylim([-20 50])
    
    subplot(4,1,2),hold on
    plot(t_sli,squeeze(d(exampl_coords(1),exampl_coords(2),:)),...
        'Color',[.3 .3 .3]);
    plot(t_sli,squeeze(d(exampl_coords(1),exampl_coords(2),:)),...
        '.','Color',[.3 .3 .3]);
    xlim([ppg_onsets(plot_nrs(1))-epoch_pre ppg_onsets(plot_nrs(end)+1)+epoch_post])
%     plot(ppg_onsets,0,'ko','MarkerSize',10)
    plot([ppg_onsets ppg_onsets],[0 .2],'k:') 
    title('T2^* signal')
%     ylim([-.8 .5])
    
    % get the upsampled (NaN) of this slice
    d_up = single(NaN(size(d,1),size(d,2),length(t_vox))); % initiate with NaNs
    d_up(:,:,ismember(round(t_vox*mux_f*srate),round(t_sli*mux_f*srate))) = d;

    % temporary response matrix for 1 slice: voxel X voxel X time X PPGonset
    temp_response_matrix = single(zeros(size(ni.data,1),size(ni.data,2),length(t),length(ppg_onsets)));
    % run through all ppg onsets
    for kk = 1:length(ppg_onsets)
        [~,ppg_find] = min(abs(t_vox-ppg_onsets(kk)));
        temp_response_matrix(:,:,:,kk) = d_up(:,:,ppg_find-floor(epoch_pre*srate_epochs):ppg_find+floor(epoch_post*srate_epochs));
    end

    t_mri = [1:size(temp_response_matrix,3)]./(mux_f*srate)-epoch_pre;
    
    % code for plotting ppg for each trial:
%     for kk = 1:length(plot_nrs)
%     plot_nr = plot_nrs(kk);
%         subplot(4,2,5),hold on
%         t_physio = [-epoch_pre*physio.ppg.srate:epoch_post*physio.ppg.srate]./physio.ppg.srate;
%         plot(t_physio,...
%             physio.ppg.data(round(ppg_onsets(plot_nr)*physio.ppg.srate)-(epoch_pre*physio.ppg.srate):...
%             round(ppg_onsets(plot_nr)*physio.ppg.srate)+(epoch_post*physio.ppg.srate)),...
%             'Color',plot_colors(kk,:))
%         xlim([-epoch_pre*physio.ppg.srate./physio.ppg.srate epoch_post*physio.ppg.srate./physio.ppg.srate])
%     end
    plot_colors = lines(8);
    
    subplot(4,2,5),hold on
    for kk = 1:length(t_mri)
        if mod(kk,2)==0
            x_fill = [t_mri(kk)-.5*mean(diff(t_mri)) t_mri(kk)+.5*mean(diff(t_mri)) t_mri(kk)+.5*mean(diff(t_mri)) t_mri(kk)-.5*mean(diff(t_mri))];
            fill(x_fill,[-1 -1 1 1],[.8 .8 .8],'EdgeColor',[.8 .8 .8])
        end
    end
    clear x_fill
    for kk = 1:length(plot_nrs) % even trials
        plot_nr = plot_nrs(kk);
        data_plot = squeeze(temp_response_matrix(exampl_coords(1),exampl_coords(2),:,plot_nr));
        plot(t_mri(~isnan(data_plot)),data_plot(~isnan(data_plot)),...
            'Color',[.3 .3 .3])
        plot(t_mri,data_plot,...
            '.','Color',[.3 .3 .3],'MarkerSize',10)%,plot_colors(kk,:))
    end
    xlim([-epoch_pre*physio.ppg.srate./physio.ppg.srate epoch_post*physio.ppg.srate./physio.ppg.srate])
%     ylim([-.8 .5])
    plot([0 0],[-.8 .5],'k:')
    
    subplot(4,2,6),hold on % odd trials
    for kk = 1:length(t_mri)
        if mod(kk,2)==0
            x_fill = [t_mri(kk)-.5*mean(diff(t_mri)) t_mri(kk)+.5*mean(diff(t_mri)) t_mri(kk)+.5*mean(diff(t_mri)) t_mri(kk)-.5*mean(diff(t_mri))];
            fill(x_fill,[-1 -1 1 1],[.8 .8 .8],'EdgeColor',[.8 .8 .8])
        end
    end
    clear x_fill
    for kk = 1:length(plot_nrs)
        plot_nr = plot_nrs(kk)+1;
        data_plot = squeeze(temp_response_matrix(exampl_coords(1),exampl_coords(2),:,plot_nr));
        plot(t_mri(~isnan(data_plot)),data_plot(~isnan(data_plot)),...
            'Color',[.3 .3 .3])
        plot(t_mri,data_plot,...
            '.','Color',[.3 .3 .3],'MarkerSize',10)%,plot_colors(kk,:))
    end
    xlim([-epoch_pre*physio.ppg.srate./physio.ppg.srate epoch_post*physio.ppg.srate./physio.ppg.srate])
%     ylim([-.8 .5])
    plot([0 0],[-.8 .5],'k:')
    
    subplot(4,2,7),hold on
    for kk = 1:length(t_mri)
        if mod(kk,2)==0
            x_fill = [t_mri(kk)-.5*mean(diff(t_mri)) t_mri(kk)+.5*mean(diff(t_mri)) t_mri(kk)+.5*mean(diff(t_mri)) t_mri(kk)-.5*mean(diff(t_mri))];
            fill(x_fill,[-1 -1 1 1],[.8 .8 .8],'EdgeColor',[.8 .8 .8])
        end
    end
    data_plot_even = nanmean(squeeze(temp_response_matrix(exampl_coords(1),exampl_coords(2),:,[2:2:end])),2);
    plot(t_mri,data_plot_even,...
        '-','Color',[0 0 0],'LineWidth',1)
    plot(t_mri,data_plot_even,...
        '.','Color',[0 0 0],'MarkerSize',10)
    xlim([-epoch_pre*physio.ppg.srate./physio.ppg.srate epoch_post*physio.ppg.srate./physio.ppg.srate])
    title('average even')
%     ylim([-.8 .5])
    plot([0 0],[-.8 .5],'k:')

    subplot(4,2,8),hold on
    for kk = 1:length(t_mri)
        if mod(kk,2)==0
            x_fill = [t_mri(kk)-.5*mean(diff(t_mri)) t_mri(kk)+.5*mean(diff(t_mri)) t_mri(kk)+.5*mean(diff(t_mri)) t_mri(kk)-.5*mean(diff(t_mri))];
            fill(x_fill,[-1 -1 1 1],[.8 .8 .8],'EdgeColor',[.8 .8 .8])
        end
    end
    data_plot_odd = nanmean(squeeze(temp_response_matrix(exampl_coords(1),exampl_coords(2),:,[1:2:end])),2);
    plot(t_mri,data_plot_odd,...
        '-','Color',[0 0 0],'LineWidth',1)
    plot(t_mri,data_plot_odd,...
        '.','Color',[0 0 0],'MarkerSize',10)
    xlim([-epoch_pre*physio.ppg.srate./physio.ppg.srate epoch_post*physio.ppg.srate./physio.ppg.srate])
    title('average odd')
%     ylim([-.8 .5])
    plot([0 0],[-.8 .5],'k:')
%     clear d_up d
    
end

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/methods/s' int2str(s_nr) '_scan' int2str(scan_nr) '_SignalTraces'])
% print('-painters','-r300','-depsc',[dDir './figures/methods/s' int2str(s_nr) '_scan' int2str(scan_nr) '_SignalTraces'])
%%
exampl_coords = [27 45 26];%[27 51 20] = superior saggital sinus, S1, scan3

sliceThisDim = 3;
curPos = mrAnatXformCoords(acpcXform,exampl_coords); % xyz coordinates
imDims = [-90 -120 -100; 90 130 110];

niFunc.data = mean(ni.data(:,:,:,1),4); % overlay the mean functional
[imgSlice,x,y,imgSlice1,x1,y1] = bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos);

figure
image(x,y,cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:))); % background
hold on
axis image
plot(curPos(1),curPos(2),'w*','MarkerSize',10)

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/methods/s' int2str(s_nr) '_scan' int2str(scan_nr) '_SignalTraces_Location'])


%%

% add to this movie by setting small R voxels to zero

% niPlot = niftiRead('/Users/m206305/Documents/data/BrainBeat/data/20141017_1242/7_1_mux8fov4_r1_25s_4mmFA48/8202_7_1_PPGtrigResponse.nii.gz');
% load /Users/m206305/Documents/data/BrainBeat/data/20141017_1242/7_1_mux8fov4_r1_25s_4mmFA48/8202_7_1_PPGtrigResponseT.mat


videoname = ['./local/CardiacPulseMovie01'];

vidObj = VideoWriter(videoname,'MPEG-4'); %

open(vidObj); 


sliceThisDim = 3;
imDims = [-90 -120 -120; 90 130 90];
curPos = [1,1,-20];
fid = figure('Position',[0 0 500 500])
for kk = 1:length(t)
    niDF = niPlot;
    niDF.data = niPlot.data(:,:,:,kk);
    bbOverlayFuncAnat(niDF,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.02,0);
    title(['t = ' num2str(t(kk))])
    set(gcf,'PaperPositionMode','auto')
    
    % Write each frame to the file.
    for m=1:3 % write X frames
        writeVideo(vidObj,getframe(fid));
    end
    
end

close(vidObj);

%%

x = physioGet(physio,'ppgppgcurve');
ppg_t = physioGet(physio,'ppgppgtcurve');

figure,plot(ppg_t,x)
plot()


