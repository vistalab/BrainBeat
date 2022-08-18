
% This script generates Figure 1 from the manuscript titled:
%
%%%%% Measuring brain beats: cardiac-aligned fast fMRI signals %%%%%
% Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell
%

clear all
close all

%% Data directory 

[~,dDir] = bbPath;


%% Load the SMS data for subject 1

sub_labels = {'1'}; 
ses_labels = {'1'}; 
acq_labels = {'4mmFA48'};
run_nrs = {[1]};

ss = 1;%:length(sub_labels) % subjects/ses/acq
rr = 1;% run_nr
sub_label = sub_labels{ss};
ses_label = ses_labels{ss};
acq_label = acq_labels{ss};

run_nr = run_nrs{ss}(rr);
        
% load the fMRI data
fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
    ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
fmri_name = fullfile(dDir,fmri_BIDSname);
ni = niftiRead(fmri_name);

save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);

% load coregistration matrix (for the functionals):
load([save_name_base '_AcpcXform_new.mat']);
acpcXform = acpcXform_new; clear acpcXform_new

% load the T1w
t1w_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'anat',...
            ['sub-' sub_label '_ses-' ses_label '_T1w.nii']);
niAnatomy = niftiRead(fullfile(dDir,t1w_BIDSname));

%% Make Figure 1: Slice with example voxel

exampl_coords = [27 51 20];

sliceThisDim = 3;
curPos = mrAnatXformCoords(acpcXform,exampl_coords); % xyz coordinates
imDims = [-90 -120 -100; 90 130 110];

% % this can show the overlay for checking data
% niFunc = ni; % initialize for plotting
% niFunc.data = mean(ni.data(:,:,:,1),4); % overlay the mean functional
% [imgSlice,x,y,imgSlice1,x1,y1] = bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos);

figure
image(x,y,cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:))); % background
hold on
axis image
plot(curPos(1),curPos(2),'w*','MarkerSize',10)

set(gcf,'PaperPositionMode','auto')    
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
    ['Figure1_sub-' int2str(sub_label) '_SignalTracesSlice']))

%% Make Figure 1: Timeseries and waveforms

plot_nrs = [13:2:20]; % heartbeat nrs to plot

figure('Position',[0 0 250 400])

slices=[1:size(ni.data,3)];

% get the nifti stuff we need:
timing = bbGet(ni,'timing');
mux_f = bbGet(ni,'super slices');
srate = 1/bbGet(ni,'tr');

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
    plot([ppg_onsets ppg_onsets],[-20 50],'k:') 
    xlim([ppg_onsets(plot_nrs(1))-epoch_pre ppg_onsets(plot_nrs(end)+1)+epoch_post])
    title('PPG signal')
    
    subplot(4,1,2),hold on
    plot(t_sli,squeeze(d(exampl_coords(1),exampl_coords(2),:)),...
        'Color',[.3 .3 .3]);
    plot(t_sli,squeeze(d(exampl_coords(1),exampl_coords(2),:)),...
        '.','Color',[.3 .3 .3]);
    xlim([ppg_onsets(plot_nrs(1))-epoch_pre ppg_onsets(plot_nrs(end)+1)+epoch_post])
    plot([ppg_onsets ppg_onsets],[0 .2],'k:') 
    title('T2^* signal')
    
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

    t_mri = (1:size(temp_response_matrix,3))./(mux_f*srate)-epoch_pre;
    
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
            '.','Color',[.3 .3 .3],'MarkerSize',10)
    end
    xlim([-epoch_pre*physio.ppg.srate./physio.ppg.srate epoch_post*physio.ppg.srate./physio.ppg.srate])
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
            '.','Color',[.3 .3 .3],'MarkerSize',10)
    end
    xlim([-epoch_pre*physio.ppg.srate./physio.ppg.srate epoch_post*physio.ppg.srate./physio.ppg.srate])
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
    plot([0 0],[-.8 .5],'k:')
    
end

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
    ['Figure1_sub-' sub_label '_ses-' ses_label '_SignalTraces']))
print('-r300','-depsc2',fullfile(dDir,'derivatives','figures',...
    ['Figure1_sub-' sub_label '_ses-' ses_label '_SignalTraces']))
