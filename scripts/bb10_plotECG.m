
sub_labels = {'4'}; 
ses_labels = {'1'}; 
acq_labels = {'4mmFA48'};
run_nrs = {[1]};

ss = 1;    
sub_label = sub_labels{ss};
ses_label = ses_labels{ss};
acq_label = acq_labels{ss};
rr = 1;    
run_nr = run_nrs{ss}(rr);
        
fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
fmri_name = fullfile(dDir,fmri_BIDSname);
ni = niftiRead(fmri_name);

% get the nifti stuff we need:
timing = bbGet(ni,'timing'); % timing per slice
mux_f = bbGet(ni,'super slices');
srate = 1/bbGet(ni,'tr');

% get physio stuff we need:
physio = physioCreate('nifti',ni);
ppg_onsets = physioGet(physio,'ppg peaks');
ecg_onsets = physioGet(physio,'ecg peaks');


%%
figure

subplot(2,1,1),hold on

ecg_srate = physio.ecg.srate;
ecg_plot = zscore(physio.ecg.data);
plot([1:length(physio.ecg.data)]/ecg_srate,ecg_plot,'Color',[1 .7 .9])
plot(ecg_onsets,1,'r*')

band = 5;
Rp   = 3; Rs = 60; % third order Butterworth
high_p =  band(1)*2/ecg_srate;
delta = 0.001*2/ecg_srate;
high_s = min(1-delta,high_p+0.1);
[n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
[bf_b,bf_a] = butter(n_band,wn_band,'low');
ecg_filt    = filtfilt(bf_b,bf_a,physio.ecg.data);
plot([1:length(physio.ecg.data)]/ecg_srate,zscore(ecg_filt),'m')

ppg_srate = physio.ppg.srate;
ppg_plot = zscore(physio.ppg.data);
plot([1:length(physio.ppg.data)]/ppg_srate,ppg_plot,'b')
plot(ppg_onsets,1,'b*')

xlim([79 92])
xlabel('time (s)')

subplot(2,1,2),hold on
hist(ppg_onsets-ecg_onsets,[.1:.005:.2],'w')
xlabel('Delta ecg and ppg onset (s)')
ylabel('number of heart beats')

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','ECGvsPPG_sub4'))
print('-painters','-r300','-depsc2',fullfile(dDir,'derivatives','brainbeat','ECGvsPPG_sub4'))


%%


signal = physio.ecg.data;
srate = physio.ecg.srate;
interval = .7*physioGet(physio,'ppg rate');

% low-pass filter 
% not necessary, keep in for lower quality data?
band = 5;
Rp   = 3; Rs = 60; % third order Butterworth
high_p =  band(1)*2/srate;
delta = 0.001*2/srate;
high_s = min(1-delta,high_p+0.1);

[n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
[bf_b,bf_a] = butter(n_band,wn_band,'low');
band_sig    = filtfilt(bf_b,bf_a,signal);


figure,hold on
ppg_srate = physio.ppg.srate;
ecg_srate = physio.ecg2.srate;
plot([1:length(physio.ppg.rawdata)]/ppg_srate,zscore(physio.ppg.rawdata))
plot([1:length(physio.ecg2.rawdata)]/ecg_srate,zscore(band_sig))
%%
% detect peak sample positions
[~,onsets] = findpeaks(band_sig,'minpeakdistance',round(interval*srate));

plot(onsets/ecg_srate,'r*');
ppg_onsets = physioGet(physio,'ppg peaks');
plot(ppg_onsets/ppg_srate,'b*');

