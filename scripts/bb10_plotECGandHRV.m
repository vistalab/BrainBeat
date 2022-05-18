
%% plot ECG

all_phys = [];

sub_labels = {'4','5','1'}; 
ses_labels = {'1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1]};

for ss = 1:length(sub_labels)    
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    rr = 1;    
    run_nr = run_nrs{ss}(rr);

    fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
    fmri_name = fullfile(dDir,fmri_BIDSname);
    ni = niftiRead(fmri_name);

    % get physio stuff we need:
    all_phys(ss).physio = physioCreate('nifti',ni);
    all_phys(ss).ppg_onsets = physioGet(all_phys(ss).physio,'ppg peaks');
    all_phys(ss).ecg_onsets = physioGet(all_phys(ss).physio,'ecg peaks');
end

% remove last ecg onset for subject 3, this is nout found in ppg, otherwise
% match
all_phys(3).ecg_onsets(end) = [];

%%
figure('Position',[0 0 800 300])

for ss = 1:length(sub_labels) 

    subplot(2,3,ss),hold on

    ecg_srate = all_phys(ss).physio.ecg.srate;
    ecg_plot = zscore(all_phys(ss).physio.ecg.data);
    plot([1:length(all_phys(ss).physio.ecg.data)]/ecg_srate,ecg_plot,'Color',[1 .7 .9])
    plot(all_phys(ss).ecg_onsets,1,'r*')

    band = 5;
    Rp   = 3; Rs = 60; % third order Butterworth
    high_p =  band(1)*2/ecg_srate;
    delta = 0.001*2/ecg_srate;
    high_s = min(1-delta,high_p+0.1);
    [n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
    [bf_b,bf_a] = butter(n_band,wn_band,'low');
    ecg_filt    = filtfilt(bf_b,bf_a,all_phys(ss).physio.ecg.data);
    plot([1:length(all_phys(ss).physio.ecg.data)]/ecg_srate,zscore(ecg_filt),'m')

    ppg_srate = all_phys(ss).physio.ppg.srate;
    ppg_plot = zscore(all_phys(ss).physio.ppg.data);
    plot([1:length(all_phys(ss).physio.ppg.data)]/ppg_srate,ppg_plot,'b')
    plot(all_phys(ss).ppg_onsets,1,'b*')

    xlim([79 92])
    xlabel('time (s)')

    subplot(2,3,3+ss),hold on
    if length(all_phys(ss).ppg_onsets)==length(all_phys(ss).ecg_onsets)
        hist(all_phys(ss).ppg_onsets-all_phys(ss).ecg_onsets,[.1:.005:.3],'w')
        xlabel('Delta ecg and ppg onset (s)')
        ylabel('number of heart beats')
    end
end
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','ECGvsPPG_sub451'))
% print('-painters','-r300','-depsc2',fullfile(dDir,'derivatives','brainbeat','ECGvsPPG_sub451'))


%% get Heart Rate Variability

all_phys = [];

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

for ss = 1:length(sub_labels)    
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    rr = 1;    
    run_nr = run_nrs{ss}(rr);

    fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
    fmri_name = fullfile(dDir,fmri_BIDSname);
    ni = niftiRead(fmri_name);

    % get physio stuff we need:
    all_phys(ss).physio = physioCreate('nifti',ni);
    all_phys(ss).ppg_onsets = physioGet(all_phys(ss).physio,'ppg peaks');
end

%% plot HRV

cl_use = lines(6);

figure

% just subsequent differences (RRs)
subplot(2,1,1),hold on
for ss = 1:6
    plot(diff(all_phys(ss).ppg_onsets),'Color',cl_use(ss,:));
end

% RMSSD
int_length = 20; % seconds to calculate RMSSD
all_hrv = zeros(6,1);
for ss = 1:6
    subplot(2,4,5:7),hold on
    these_int = 1:int_length:max(all_phys(ss).ppg_onsets)-int_length;
    this_hrv = zeros(length(these_int),1);
    for kk = 1:length(these_int)
        % get 10 seconds
        these_beats = all_phys(ss).ppg_onsets(all_phys(ss).ppg_onsets>these_int(kk) & all_phys(ss).ppg_onsets<=these_int(kk)+int_length);
        % get RR (peak to peak difference),
        RR = diff(these_beats);
        % difference in ms between RRs
        these_diff = diff(RR)*1000; 
        this_hrv(kk) = sqrt(mean(these_diff.^2));
    end
    plot(this_hrv,'Color',cl_use(ss,:))
    all_hrv(ss) = mean(this_hrv);
    
    subplot(2,4,8),hold on
    b = bar(ss,all_hrv(ss));
    set(b,'FaceColor',cl_use(ss,:))
end


% b.CData = cl_use;

