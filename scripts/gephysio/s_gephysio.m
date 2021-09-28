%% s_gephysio
%
% Author: Hua Wu
%
% Illustrates how to download data from the iris project, part of the
% LEANEW1 group 
% 
% Then we read the data.
% Then we process the data using the DH methods.
% Then we save the file.
%
% We will convert this into a Gear at some point.
%

% General information so we can find the physio data
group         = 'leanew1';
project_label = 'iris';
subj_label    = 'IRIS032';
session_label = '23792';
acq_label_rest = 'resting_2';
acq_label_fmap = 'spiral high-res fieldmap';
acq_label_anat = 'T1w .8mm sag';

group_StanfordLabs         = 'cni';
project_label_StanfordLabs = 'FreeSurferTest';
subj_label_StanfordLabs    = 'IRIS032';
session_label_StanfordLabs = 'brainbeats';
ana_label_StanfordLabs     = 'freesurfer-recon-all';

fn_anat = 't1w';
fn_rest = 'rest';
fn_fmap_mag = 'fmap_mag';
fn_fmap_frq = 'fmap_frq';
fn_FS_anat = 'FS_T1';
fn_FS_aseg = 'FS_aparcDKT';

%% Here are the scitran commands for finding and downloading the data

% see bbResponse2physio in BrainBeats

% Let's try with scitran for a while.
% st = scitran('cni');

% Or use Flywheel MATLAB SDK. fw = flywheel.Flywheel(apikey)
if exist('fw', 'var') && isa(fw, 'flywheel.Flywheel')
    userInfo = fw.getCurrentUser;
    if ~isempty(userInfo)
        fprintf('Continue with current Flywheel user: %s \n', userInfo.id);
    else
        warning('Empty Flywheel user info. \n');
    end
else
    fprintf('Connect to Flywheel... \n');
    apikey = input('Please enter apikey: ', 's');
    fw = flywheel.Flywheel(apikey);
    userInfo = fw.getCurrentUser;
    fprintf('Logged into Flywheel as %s\n', userInfo.id);
end
localdir = fullfile(bbPath,'local');
datadir = fullfile(localdir,subj_label);
if ~exist(datadir, 'dir')
    mkdir(datadir);
end

%% Download files from Flywheel using scitran - resting state and physio

% Using lookup to find data. Examples of lookup methods that include the subject
lookup = fullfile(group,project_label,subj_label,session_label,acq_label_rest);
thisAcq = fw.lookup(lookup);

% See the file names
for f=1:numel(thisAcq.files) 
    % fprintf([thisAcq.files{f}.name '\n']); 
    if contains(thisAcq.files{f}.name, 'gephysio')
        gephysioZip = thisAcq.files{f};
    elseif contains(thisAcq.files{f}.name, 'dicom')
        dcmFile = thisAcq.files{f};
    elseif contains(thisAcq.files{f}.name, 'nii')
        niFile = thisAcq.files{f};
    end
end


% Download and unzip the physio files
gephysioFile = fullfile(datadir, gephysioZip.name);
gephysioZip.download(gephysioFile);
physiofiles = unzip(gephysioFile, datadir);
niFile.download(fullfile(datadir, niFile.name));

% set physio data filenames and sampling rates
param.ppg.dt  = 10;       % PPG data samples at 10ms
param.resp.dt = 40;       % Respiratory data samples  at 40 ms
param.ppg.wave.fn  = physiofiles{contains(physiofiles,'PPGData')};
param.ppg.trig.fn  = physiofiles{contains(physiofiles,'PPGTrig')};
param.resp.wave.fn = physiofiles{contains(physiofiles,'RESPData')};
param.resp.trig.fn = physiofiles{contains(physiofiles,'RESPTrig')};
% Get the acquisition parameters from the dicom, nifti file metadata (info)
param.TR             = dcmFile.info.RepetitionTime;
param.nvols          = dcmFile.info.NumberOfTemporalPositions;
param.scan_duration  = param.TR * param.nvols;   % milliseconds
param.sliceTiming    = niFile.info.SliceTiming * 1000;  % converting to milliseconds
[~,param.sliceOrder] = sort(param.sliceTiming);
param.echoSpacing    = niFile.info.EffectiveEchoSpacing;

%% download anatomical and fieldmap
anatAcq = fw.lookup(fullfile(group,project_label,subj_label,session_label,acq_label_anat));

for f=1:numel(anatAcq.files) 
    if contains(anatAcq.files{f}.name, '.nii.gz')
        anatNiFile = anatAcq.files{f};
    end
end
anatNiFile.download(fullfile(datadir, anatNiFile.name));

fmapAcq = fw.lookup(fullfile(group,project_label,subj_label,session_label,acq_label_fmap));
for f=1:numel(fmapAcq.files) 
    if contains(fmapAcq.files{f}.name, '_1.nii.gz')
        fmap_mag = fmapAcq.files{f};
    elseif contains(fmapAcq.files{f}.name, '_fieldmap.nii.gz')
        fmap_frq = fmapAcq.files{f};
    end
end
fmap_mag.download(fullfile(datadir, fmap_mag.name));
fmap_frq.download(fullfile(datadir, fmap_frq.name));

%% prepare the nifti -- reorient to standard orientation
system(['/usr/local/fsl/bin/fslreorient2std ' fullfile(datadir,niFile.name) ' ' fullfile(datadir,fn_rest)]);
system(['/usr/local/fsl/bin/fslreorient2std ' fullfile(datadir,fmap_mag.name) ' ' fullfile(datadir,fn_fmap_mag)]);
system(['/usr/local/fsl/bin/fslreorient2std ' fullfile(datadir,fmap_frq.name) ' ' fullfile(datadir,fn_fmap_frq)]);
system(['/usr/local/fsl/bin/fslreorient2std ' fullfile(datadir,anatNiFile.name) ' ' fullfile(datadir,fn_anat)]);


%% Download Freesurfer analysis results from Flywheel stanfordlabs site
if exist('fw_StanfordLabs', 'var') && isa(fw_StanfordLabs, 'flywheel.Flywheel')
    userInfo = fw_StanfordLabs.getCurrentUser;
    if ~isempty(userInfo)
        fprintf('Continue with current Flywheel user: %s \n', userInfo.id);
    else
        warning('Empty Flywheel user info. \n');
    end
else
    fprintf('Connect to Flywheel... \n');
    apikey = input('Please enter apikey: ', 's');
    fw_StanfordLabs = flywheel.Flywheel(apikey);
    userInfo = fw_StanfordLabs.getCurrentUser;
    fprintf('Logged into Flywheel as %s\n', userInfo.id);
end

ana=fw_StanfordLabs.lookup(fullfile(group_StanfordLabs, project_label_StanfordLabs, subj_label_StanfordLabs, session_label_StanfordLabs, ...
    'analyses', ana_label_StanfordLabs));
for f=1:numel(ana.files) 
    if strcmp(ana.files{f}.name(end-3:end),'.zip'), break; end
end
archive=ana.files{f};
destFile='freesurfer.zip';
archive.download(fullfile(datadir,destFile));
unzip(fullfile(datadir,destFile),datadir);

% archive.downloadZipMember('T1.mgz', [subj_label, '/mri/T1.mgz']);
% archive.downloadZipMember('aparcDKT.mgz', [subj_label, '/mri/aparc.DKTatlas+aseg.mgz']);

system(['export FREESURFER_HOME=/Applications/freesurfer; /Applications/freesurfer/bin/mri_convert ' ...
    fullfile(datadir,subj_label,'mri/T1.mgz') ' ' fullfile(datadir,[fn_FS_anat '.nii.gz'])]);
system(['export FREESURFER_HOME=/Applications/freesurfer; /Applications/freesurfer/bin/mri_convert ' ...
    fullfile(datadir,subj_label,'mri/aparc.DKTatlas+aseg.mgz') ' ' fullfile(datadir,[fn_FS_aseg '.nii.gz'])]);


%% Read the physio data and find peaks

physioType = {'ppg', 'resp'};

% First read in the waveform data to determine the total length of the recording,
% then trim the trigger data accordingly
for p = 1:numel(physioType)
    [param.(physioType{p}).wave.data_raw, param.(physioType{p}).wave.data_sync, ...
        param.(physioType{p}).wave.t_raw, param.(physioType{p}).wave.t_sync, param.scanStart] = ...
        physioRead(param.(physioType{p}).wave.fn, param.(physioType{p}).dt, param.scan_duration, 'wave');
end

for p = 1:numel(physioType)
    [param.(physioType{p}).trig.data_raw, param.(physioType{p}).trig.data_sync,~,~,~] = ...
        physioRead(param.(physioType{p}).trig.fn, param.(physioType{p}).dt, param.scan_duration, 'trig', param.scanStart);
end

% Filter signal and find peaks in ppg and resp data using Dora's function
for p = 1:numel(physioType)
    [param.(physioType{p}).trig.data_filt, param.(physioType{p}).wave.data_filt] = ...
        physioPeaks(param.(physioType{p}).wave.data_sync, param.(physioType{p}).dt);
end

%% Plot a segment of the data

plot_window = 10000;    % Duration of physio data to plot for inspection (10 seconds)
plot_start  = 4000; 
if plot_start<param.resp.dt, plot_start=param.resp.dt; end
if plot_start>param.scan_duration-plot_window, plot_start=param.scan_duration-plot_window; end
plot_end    = plot_start + plot_window;

figure; 
for p = 1:numel(physioType)
    subplot(2,1,p);
    r = (plot_start < param.(physioType{p}).wave.t_sync) & (param.(physioType{p}).wave.t_sync < plot_end);
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_sync(r),'b'); hold on;
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_filt(r),'r:', 'LineWidth',1); 
    
    trig_in_window = param.(physioType{p}).trig.data_filt(...
        (plot_start < param.(physioType{p}).trig.data_filt) & ...
        (param.(physioType{p}).trig.data_filt < plot_end));
    plot(trig_in_window, ones(size(trig_in_window))*max(param.(physioType{p}).wave.data_sync),'ro');
    getrig_in_window = param.(physioType{p}).trig.data_sync(...
        (plot_start < param.(physioType{p}).trig.data_sync) & ...
        (param.(physioType{p}).trig.data_sync < plot_end));
    plot(getrig_in_window, ones(size(getrig_in_window))*max(param.(physioType{p}).wave.data_sync),'bx');
    xlabel('Time (ms)'); ylabel(sprintf('%s recorded data', physioType{p})); xlim([plot_start, plot_end]);
end

%% find ppg trigger for each fmri volume
% 
% param.rsfmri.t = (0:param.TR:(param.scan_duration - param.TR));
% hbeats = param.ppg.trig.data_filt;
% for sl = 1:size(rsfmri, 3)
%     for ii = 1:numel(param.rsfmri.t)
%         tt = param.rsfmri.t(ii) + param.sliceTiming(sl);    % correct for sliceTiming for each slice
%         hbeats_before = (hbeats < tt);
%         this_beat = hbeats(find(hbeats_before, 1, 'last'));
%         if isempty(this_beat)
%             param.rsfmri.beat_tdiff(sl,ii) = NaN;
%         else
%             param.rsfmri.beat_tdiff(sl,ii) = tt - this_beat;
%         end
%     end
% end

%% preprocessing the fmri data
% fieldmap correction & alignment of the brain atlas
cd(datadir);
system(['export FSLDIR=/usr/local/fsl;' ...
    fullfile(bbPath,'scripts/fsl-commands.sh') ' ' fn_anat ' ' fn_rest ' ' fn_FS_anat ' ' fn_FS_aseg ' ' fn_fmap_mag ' ' fn_fmap_frq ' ' sprintf('%.6f',param.echoSpacing) ]);
cd(bbPath);

% detrend the fmri time series 
rsfmri = niftiread(fullfile(datadir,[fn_rest '_unwarp.nii.gz']));
dummy_cycles = 3;

rsfmri_reg = physioDetrend(rsfmri, dummy_cycles);

%% find the fmri volume closest to each ppg trigger

[y,x,z] = meshgrid(1:size(rsfmri_reg,1),1:size(rsfmri,2),1:size(rsfmri,3));  % x,y coordinates are swapped in meshgrid
roi = [reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
hbeats = param.ppg.trig.data_filt;

[gatedSignal,gatedTime,gatedSignalNormalized,gatedTimesNormalized] = physioResponse(hbeats, rsfmri_reg, 'roi',roi);

%% average the ppg response over a certain time window
timeAvgWindow = 50;  % in ms
timeResolution = 10;
[avgGatedSignal,interpTimes] = physioResponseAvg(gatedSignal, gatedTime, 'timeAvgWindow',timeAvgWindow,'timeResolution',timeResolution);
[avgGatedSignalNormalized,interpTimesNormalized] = physioResponseAvg(gatedSignalNormalized, gatedTimesNormalized, 'timeAvgWindow',timeAvgWindow,'timeResolution',timeResolution);


%% compute average response in ROIs defined in DKT atlas
dkt_table = readtable(fullfile(bbPath,'colormaps/dkt_areas.tsv'),'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dkt_table_surface = readtable(fullfile(bbPath,'colormaps/dkt_areas_surface.tsv'),'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dkt_atlas = niftiread(fullfile(datadir,'fs_aparc_resampled.nii.gz'));
roiNames = dkt_table.label;
roiCodes = dkt_table.label_nr;

% reshape DKT atlas to a vector across the spatial dimensions
dkt_atlas_v = reshape(dkt_atlas,[],1);
avgROISig = zeros(length(roiNames), size(avgGatedSignal,2));
avgROISigNormalized = zeros(length(roiNames), size(avgGatedSignalNormalized,2));
rsfmri_v = reshape(rsfmri,[],size(rsfmri,4));
% compute the average in each ROI
for kk = 1:length(roiNames)
    avgROIraw(kk) = mean(rsfmri_v(dkt_atlas_v==roiCodes(kk),dummy_cycles+1),1);
    avgROISig(kk,:) = mean(avgGatedSignal(dkt_atlas_v==roiCodes(kk),:),1) / avgROIraw(kk) * 100; 
    avgROISigNormalized(kk,:) = mean(avgGatedSignalNormalized(dkt_atlas_v==roiCodes(kk),:),1) / avgROIraw(kk) * 100; 
end

%% look at some ROIs with high signal variation
mean_hrate = mean(diff(hbeats))/1000;
for kk = 1:size(avgROISig,1)
    % variation = max(avgROISig(kk,:))-min(avgROISig(kk,:)); 
    variation = max(avgROISig(kk,:))-min(avgROISig(kk,:)); 
    if variation > 1
        fprintf('ROI %d %s, value %d\n',kk, roiNames{kk}, roiCodes(kk)); 
        figure; hold on;
        plot(interpTimes,avgROISig(kk,:),'DisplayName',sprintf('%s', roiNames{kk}));
        plot(interpTimesNormalized*mean_hrate,avgROISigNormalized(kk,:),'DisplayName',sprintf('%s normalized', roiNames{kk}));
        xlim([0,param.TR]); xlabel('Time (ms)'); ylabel('Contrast (%)');
        ax=gca; ax.FontSize=16; legend show; legend('FontSize',16); 
    end
end

%% plot a specific ROI by name
% roiname = 'rightentorhinal';
% for kk = 1:length(roiNames)
%     if strcmp(roiNames{kk}, roiname)
%         figure; hold on;
%         plot(interpTimes,avgROISigNormalized(kk,:),'DisplayName',sprintf('%s', roiNames{kk}));
%         xlim([0,2000]); xlabel('Time (ms)'); ylabel('Contrast (%)');
%         ax=gca; ax.FontSize=16; legend show; legend('FontSize',16); 
%     end
% end

%% SMS BOLD data
% load('/Users/huawu/Documents/GitHub/BrainBeat/local/sms8/sub-1_ses-1_acq-4mmFA48_run-1_PPGtrigResponseT.mat', 't')
% sms_response=niftiread('local/sms8/sub-1_ses-1_acq-4mmFA48_run-1_PPGtrigResponse.nii.gz');
% size(sms_response)
% sms_atlas=niftiread('local/sms8/sub-1_ses-1_acq-4mmFA48_run-1_r_DKTatlas_aseg.nii.gz');
% sms_response_v=reshape(sms_response,[],size(sms_response,4));
% sms_atlas_v=reshape(sms_atlas,[],1);
% for kk=1:length(roiNames),sms_avgResponse(kk,:)=mean(sms_response_v(sms_atlas_v==roiCodes(kk),:),1);end
% figure;hold on; 
% for kk=1:length(roiNames)
%     v=max(sms_avgResponse(kk,:))-min(sms_avgResponse(kk,:)); 
%     if v>0.02 
%         plot(sms_t,sms_avgResponse(kk,:),'DisplayName',sprintf('%s',roiNames{kk}));
%     end
% end;legend show
% 
% kk=61;
% figure;hold on; 
% plot(sms_t,sms_avgResponse(kk,:),'DisplayName','sms');
% plot(sms_t/1.25*0.78,sms_avgResponse(kk,:),'k--','DisplayName','sms scaled');
% plot(((1:size(avgROISig,2))+2)*0.01,avgROISig(kk,:)/avgROIraw(kk),'r','DisplayName','non-sms');
% title(roiNames{kk},'FontSize',16);legend('FontSize',16);xlabel('t (s)','FontSize',16);ylabel('contrast','FontSize',16)

%% Save output
% reshapde the outputs to 4D -- only keep data points within 2.5s
avgGatedSignalNormalizedC = avgGatedSignalNormalized(:,interpTimesNormalized<2500);
interpTimesNormalizedC = interpTimesNormalized(interpTimesNormalized<2500);
avgGatedSignalNormalized4D = reshape(avgGatedSignalNormalizedC,[size(rsfmri,1:3),size(avgGatedSignalNormalizedC,2)]);
% save avgGatedSignal into 4D NIFTI
ni_info = niftiinfo(fullfile(datadir,fn_rest)) ;
ni_info.Datatype = 'double';
ni_info.ImageSize = size(avgGatedSignalNormalized4D);
niftiwrite(avgGatedSignalNormalized4D, fullfile(datadir,'avgGatedSignalNormalized'), ni_info, 'Compressed', true);
save(fullfile(datadir,'avgTimes.mat'),'interpTimesNormalizedC');
% % average heart beat / respiration rate in Hz
% for p = 1:numel(physioType)
%     dtrig = diff(param.(physioType{p}).trig.data_filt);
%     param.(physioType{p}).trig.average_rate = 1000 / mean(dtrig(abs(dtrig - mean(dtrig)) < 3 * std(dtrig))); 
% end
% 
% fid = fopen(fullfile(datadir, 'physio.json'), 'w');
% fprintf(fid, jsonencode(param, 'PrettyPrint', true));
% fclose(fid);

%% END