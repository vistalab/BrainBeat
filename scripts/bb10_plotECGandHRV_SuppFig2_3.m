
% This script generates Supplemental Figures X and Y from the manuscript titled:
%
% Measuring brain beats: cardiac-aligned fast fMRI signals
% Dora Hermes, Hua Wu, Adam B. Kerr, Brian Wandell

clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%% Supplemental Figure 2
%%
%% Load ECG and PPG data to compare ECG vs PPG onsets for sub-4_ses-1 and sub-1_ses-2
%% 

all_phys = [];

sub_labels = {'4','1'};%{'4','5','1'}; 
ses_labels = {'1','2'};%{'1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48'};
run_nrs = {1,1};
% we also collected ECG data in subject 5, session 1, run 1, but that was
% extremely noisy from scanner artifact...

for ss = 1:length(sub_labels)    
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    rr = 1;    
    run_nr = run_nrs{ss}(rr);

    fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
    fmri_name = fullfile(dDir,fmri_BIDSname);
    ni = niftiRead(fmri_name);

    % get physio stuff we need:
    all_phys(ss).physio = physioCreate('nifti',ni);
    all_phys(ss).ppg_onsets = physioGet(all_phys(ss).physio,'ppg peaks');
    all_phys(ss).ecg_onsets = physioGet(all_phys(ss).physio,'ecg peaks');
end

% remove last ecg onset for subject 1, session 2, this is not found in ppg, otherwise
% all onsets match
all_phys(find(ismember(sub_labels,'1'))).ecg_onsets(end) = [];

%% Supplemental Figure 2: plot ECG and PPG signals

figure('Position',[0 0 500 300])

for ss = 1:length(sub_labels) 

    subplot(2,2,ss),hold on

    ecg_srate = all_phys(ss).physio.ecg.srate;
    ecg_plot = zscore(all_phys(ss).physio.ecg.data);
    plot((1:length(all_phys(ss).physio.ecg.data))/ecg_srate,ecg_plot,'Color',[1 .7 .9])

    band = 5;
    Rp   = 3; Rs = 60; % third order Butterworth
    high_p =  band(1)*2/ecg_srate;
    delta = 0.001*2/ecg_srate;
    high_s = min(1-delta,high_p+0.1);
    [n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
    [bf_b,bf_a] = butter(n_band,wn_band,'low');
    ecg_filt    = filtfilt(bf_b,bf_a,all_phys(ss).physio.ecg.data);
    ecg_plot2 = zscore(ecg_filt);
    plot((1:length(all_phys(ss).physio.ecg.data))/ecg_srate,ecg_plot2,'r')
    plot(all_phys(ss).ecg_onsets,ecg_plot2(round(all_phys(ss).ecg_onsets*ecg_srate)),'r*')

    ppg_srate = all_phys(ss).physio.ppg.srate;
    ppg_plot = zscore(all_phys(ss).physio.ppg.data);
    plot((1:length(all_phys(ss).physio.ppg.data))/ppg_srate,ppg_plot,'b')
    plot(all_phys(ss).ppg_onsets,ppg_plot(round(all_phys(ss).ppg_onsets*ppg_srate)),'b*')

    xlim([82 90])
    xlabel('time (s)')

    subplot(2,2,2+ss),hold on
    if length(all_phys(ss).ppg_onsets)==length(all_phys(ss).ecg_onsets)
        hist(all_phys(ss).ppg_onsets-all_phys(ss).ecg_onsets,[.1:.005:.3],'w')
        xlabel('$${\Delta}$$ ppg-ecg onsets (s)','Interpreter','latex')

        ylabel('Number of heart beats','Interpreter','latex')
        
        delta_onset = (all_phys(ss).ppg_onsets-all_phys(ss).ecg_onsets); 
        pd = fitdist(delta_onset,'Normal'); % fit a normal distribution
        x_values = 0.050:0.001:0.350;
        y = pdf(pd,x_values);
        plot(x_values,y,'LineWidth',2,'Color',[0 .5 1])
        ci95 = paramci(pd,'Alpha',.05); % returns ci of the mean and std
        title(['mean 95 \% ci ' num2str(ci95(:,1)',3)],'Interpreter','latex')
    end
end

set(gcf,'PaperPositionMode','auto')    
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures','SuppFigureS2_ECGvsPPG'))
print('-r300','-depsc2',fullfile(dDir,'derivatives','figures','SuppFigureS2_ECGvsPPG'))


%% Supplemental Figure 3
%%
%% Calculate Heart Rate Variability (HRV) for all subjects
%%

all_phys = [];

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,1,1,1,1};

for ss = 1:length(sub_labels)    
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    rr = 1;    
    run_nr = run_nrs{ss}(rr);

    fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
    fmri_name = fullfile(dDir,fmri_BIDSname);
    ni = niftiRead(fmri_name);

    % get physio stuff we need:
    all_phys(ss).physio = physioCreate('nifti',ni);
    all_phys(ss).ppg_onsets = physioGet(all_phys(ss).physio,'ppg peaks');
end

%% Supplemental Figure 3: plot HRV
cl_use = lines(6); % colors for different subjects

figure('Position',[0 0 550 300])

% just subsequent differences (RRs)
for ss = 1:6 % subjects
    subplot(2,4,1:3),hold on
    plot(all_phys(ss).ppg_onsets(2:end),diff(all_phys(ss).ppg_onsets),'Color',cl_use(ss,:));
    
    subplot(2,4,4),hold on    
    b = bar(ss,mean(diff(all_phys(ss).ppg_onsets)));
    set(b,'FaceColor',cl_use(ss,:))
    eb = errorbar(ss,mean(diff(all_phys(ss).ppg_onsets)),std(diff(all_phys(ss).ppg_onsets)));
    set(eb,'Color','k')   
end

subplot(2,4,1:3),hold on
ylabel('P-P interval (s)')
xlabel('Time (s)')
xlim([0 max(all_phys(1).ppg_onsets)]),ylim([0 1.5])

subplot(2,4,4),hold on
set(gca,'XTick',1:6),ylim([0 1.5])
xlabel('Subject #'),ylabel('P-P interval (s)')

% Root Mean Sum of Squared Differences (RMSSD)
int_length = 20; % interval in seconds to calculate RMSSD
all_hrv = zeros(6,1);
for ss = 1:6 % subjects
    subplot(2,4,5:7),hold on
    these_int = 1:int_length:max(all_phys(ss).ppg_onsets);
    this_hrv = zeros(length(these_int),1);
    for kk = 1:length(these_int)
        % get 10 seconds
        these_beats = all_phys(ss).ppg_onsets(all_phys(ss).ppg_onsets>these_int(kk) & all_phys(ss).ppg_onsets<=these_int(kk)+int_length);
        % get RR (peak to peak difference)
        RR = diff(these_beats);
        % difference in ms between RRs
        these_diff = diff(RR)*1000; 
        this_hrv(kk) = sqrt(mean(these_diff.^2));
    end
    plot(these_int-1+int_length,this_hrv,'Color',cl_use(ss,:))
    all_hrv(ss) = mean(this_hrv);
    
    subplot(2,4,8),hold on
    b = bar(ss,all_hrv(ss));
    set(b,'FaceColor',cl_use(ss,:))   
    eb = errorbar(ss,all_hrv(ss),std(this_hrv));
    set(eb,'Color','k')   
end

subplot(2,4,5:7),hold on
title('Root Mean Square of Successive Differences between heartbeats')
xlabel('Time (s)'),ylabel('RMSSD (ms)')
xlim([0 max(all_phys(1).ppg_onsets)]),ylim([0 100])
subplot(2,4,8),hold on
set(gca,'XTick',1:6),ylim([0 100])
ylabel('RMSSD (ms)')
xlabel('Subject #')

print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures','SuppFigureS3_HRV'))
print('-r300','-depsc2',fullfile(dDir,'derivatives','figures','SuppFigureS3_HRV'))
