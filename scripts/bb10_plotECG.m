
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


%%
figure,hold on
ppg_srate = physio.ppg.srate;
ecg_srate = physio.ecg2.srate;
plot([1:length(physio.ppg.rawdata)]/ppg_srate,zscore(physio.ppg.rawdata))
plot([1:length(physio.ecg2.rawdata)]/ecg_srate,zscore(physio.ecg2.rawdata-physio.ecg3.rawdata))

xlim([9 11.5])
% a = physio.ecg3.rawdata-physio.ecg2.rawdata;
% [band_sig] = ieeg_butterpass(a,[1 10],1000);
% plot([1:length(physio.ecg2.rawdata)]/ecg_srate,zscore(band_sig))

% plot([1:length(physio.ecg2.rawdata)]/ecg_srate,zscore(physio.ecg2.rawdata-physio.ecg3.rawdata))

% plot([1:length(physio_output.ecg2.rawdata)]/ecg_srate,zscore(physio_output.ecg2.rawdata))
% plot([1:length(physio_output.ecg3.rawdata)]/ecg_srate,zscore(physio_output.ecg3.rawdata))

%%
