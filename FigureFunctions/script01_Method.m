clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';
% chdir(dDir)

%% The T2* data are here.  

% Select a subject and scan nummer
s_nr = 2;
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

%% Quick overlay between functionals and anatomy

sliceThisDim = 3;
if s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [1,10,-20];
%     curPos = [-29,-14,-50]; % 03
%     curPos = [-3,30,-43]; % 03
    curPos = [1,10,-20];
    curPos = [1,1,-20];
    curPos = [-11 34 -71]; % Carotid
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
%     curPos = [0,4,38];
elseif s_nr == 4
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
end
niFunc = ni;

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

exampl_coords = [27 51 20];

figure

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
epoch_post = 2;%sec post-onset

step_size = 1/srate/mux_f;% in s
srate_epochs = 1/step_size;
t = [-epoch_pre:step_size:epoch_post];%s timing for 1 epoch
t_vox = min(timing(:)):step_size:max(timing(:)); % maximal timing accuracy for this scan session

% only use those ppg onsets that have an entire trial before and after
ppg_onsets = ppg_onsets((ppg_onsets-epoch_pre)>1); % get rid of early ones, and the first scans
ppg_onsets = ppg_onsets((ppg_onsets+epoch_post)<max(timing(:))); % get rid of late ones

% initiate output matrix to fill:
response_matrix = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));
response_matrix_std = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));
response_matrix_odd = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));
response_matrix_even = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));

figure

for s = exampl_coords(3)%1:length(slices)
    disp(['slice ' int2str(s) ' of ' int2str(length(slices))])
    sli = slices(s);
    d = squeeze(ni.data(:,:,sli,1:end));

    % demean and express in percent modulation 
    d_norm=reshape(d,[size(d,1) * size(d,2), size(d,3)]);
    points_use=4:size(d_norm,2); % do not use the first couple of scans
    for k=1:size(d_norm,1)
        % do not use first scans:
        x = points_use;
        y = d_norm(k,points_use);
        
        % detrend
        p = polyfit(x,y,1);    

        % percent modulation
        mean_factor = mean(d_norm(k,points_use)); % mean
        d_norm(k,:) = (d_norm(k,:) - (p(1)*[1:size(d_norm,2)] + p(2)))./mean_factor;
    end
    d = reshape(d_norm,[size(d,1), size(d,2), size(d,3)]);
    clear d_norm

    % get the timing for this slice
    t_sli = timing(sli,:);

    subplot(4,1,1),hold on
    plot([1:length(physio.ppg.data)]./physio.ppg.srate,physio.ppg.data,'k')
    plot(ppg_onsets,physio.ppg.data(round(ppg_onsets*physio.ppg.srate)),'ro')

    subplot(4,1,2),hold on
    plot(t_sli,squeeze(ni.data(exampl_coords(1),exampl_coords(2),exampl_coords(3),:)),'k');

    
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
    
    plot_colors = lines(8);
    plot_nrs = [51:58];
    for kk = 1:length(plot_nrs)
        plot_nr = plot_nrs(kk);
        subplot(4,1,3),hold on
        t_physio = [-epoch_pre*physio.ppg.srate:epoch_post*physio.ppg.srate]./physio.ppg.srate;
        plot(t_physio,...
            physio.ppg.data(round(ppg_onsets(plot_nr)*physio.ppg.srate)-(epoch_pre*physio.ppg.srate):...
            round(ppg_onsets(plot_nr)*physio.ppg.srate)+(epoch_post*physio.ppg.srate)),...
            'Color',plot_colors(kk,:))
        xlim([-epoch_pre*physio.ppg.srate./physio.ppg.srate epoch_post*physio.ppg.srate./physio.ppg.srate])

        subplot(4,1,4),hold on
        t_mri = [1:size(temp_response_matrix,3)]./(mux_f*srate)-epoch_pre;
        plot(t_mri,...
            squeeze(temp_response_matrix(exampl_coords(1),exampl_coords(2),:,plot_nr)),...
            '.','Color',plot_colors(kk,:))
        xlim([-epoch_pre*physio.ppg.srate./physio.ppg.srate epoch_post*physio.ppg.srate./physio.ppg.srate])
    end
    
    clear d_up d
    
    % output:
    response_matrix(:,:,s,:) = nanmean(temp_response_matrix,4);
    response_matrix_std(:,:,s,:) = nanstd(temp_response_matrix,[],4);
    response_matrix_odd(:,:,s,:) = nanmean(temp_response_matrix(:,:,:,1:2:end),4);
    response_matrix_even(:,:,s,:) = nanmean(temp_response_matrix(:,:,:,2:2:end),4);
    
    clear temp_response_matrix
end
