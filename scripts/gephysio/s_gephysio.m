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
subj_label    = 'IRIS021';
session_label = '22776';
acq_label     = 'resting_2';

%% Here are the scitran commands for finding and downloading the data

% see bbResponse2physio in BrainBeats

% Let's try with scitran for a while.
% st = scitran('cni');

% Or use Flywheel MATLAB SDK. fw = flywheel.Flywheel(apikey)
javaaddpath('~/MATLAB/flywheel-sdk/api/rest-client.jar');
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

%% Download files from Flywheel using scitran

% Using lookup to find data. Examples of lookup methods that include the subject
lookup = fullfile(group,project_label,subj_label,session_label,acq_label);
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

nii = strsplit(niFile.name, '.');
datadir = fullfile(localdir,nii{1});
if ~exist(datadir, 'dir')
    mkdir(datadir);
end

% Download and unzip the physio files
gephysioFile = fullfile(datadir, gephysioZip.name);
gephysioZip.download(gephysioFile);
physiofiles = unzip(gephysioFile, datadir);
niFile.download(fullfile(datadir, niFile.name));
rsfmri = niftiread(fullfile(datadir, niFile.name));
rsfmri_info = niftiinfo(fullfile(datadir, niFile.name));

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

param.rsfmri.t = (0:param.TR:(param.scan_duration - param.TR));
hbeats = param.ppg.trig.data_filt;
for sl = 1:size(rsfmri, 3)
    for ii = 1:numel(param.rsfmri.t)
        tt = param.rsfmri.t(ii) + param.sliceTiming(sl);    % correct for sliceTiming for each slice
        hbeats_before = (hbeats < tt);
        this_beat = hbeats(find(hbeats_before, 1, 'last'));
        if isempty(this_beat)
            param.rsfmri.beat_tdiff(sl,ii) = NaN;
        else
            param.rsfmri.beat_tdiff(sl,ii) = tt - this_beat;
        end
    end
end

%% regress the fmri time series 
dummy_cycles = 3;

rsfmri_reg = physioDetrend(rsfmri, dummy_cycles);

%% find the fmri volume closest to each ppg trigger

[y,x,z] = meshgrid(1:size(rsfmri_reg,1),1:size(rsfmri,2),1:size(rsfmri,3));  % x,y coordinates are swapped in meshgrid
roi = [reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
hbeats = param.ppg.trig.data_filt;

[gatedSignal,gatedTime] = physioResponse(hbeats, rsfmri_reg, 'roi',roi);

%% average the ppg response over a certain time window
timeWindow = 50;  % in ms
timeResolution = 10;
[avgGatedSignal,interpTimes] = physioResponseAvg(gatedSignal, gatedTime, 'timeWindow',timeWindow,'timeResolution',timeResolution);


%% plot signal from an example voxel
% x1=60; y1=39; z1=20;
% idx = sub2ind(size(rsfmri,1:3), x1,y1,z1);
% figure;  plot(interpTimes, squeeze(avgGatedSignal(idx,:)),'-x');
% hold on; plot(squeeze(gatedTime(idx,:)), squeeze(gatedSignal(idx,:)),'x');

%% compute average response in ROIs defined in DKT atlas
dkt_table = readtable(fullfile(bbPath,'colormaps/dkt_areas.tsv'),'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dkt_table_surface = readtable(fullfile(bbPath,'colormaps/dkt_areas_surface.tsv'),'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dkt_atlas = niftiread(fullfile(datadir,'aparc.DKTatlas+aseg_resampled.nii.gz'));
roiNames = dkt_table.label;
roiCodes = dkt_table.label_nr;

% reshape DKT atlas to a vector across the spatial dimensions
dkt_atlas_v = reshape(dkt_atlas,[],1);
avgROISig = zeros(length(roiNames), size(avgGatedSignal,2));
rsfmri_v = reshape(rsfmri,[],size(rsfmri,4));
% compute the average in each ROI
for kk = 1:length(roiNames)
    avgROISig(kk,:) = mean(avgGatedSignal(dkt_atlas_v==roiCodes(kk),:),1); 
    avgROIraw(kk) = mean(rsfmri_v(dkt_atlas_v==roiCodes(kk),dummy_cycles+1),1);
end

%% look at some ROIs with high signal variation
figure; hold on;
for kk = 1:size(avgROISig,1)
    % variation = max(avgROISig(kk,:))-min(avgROISig(kk,:)); 
    variation = (max(avgROISig(kk,:))-min(avgROISig(kk,:)))/avgROIraw(kk); 
    if variation > 0.010 && variation < 0.011
        fprintf('ROI %d %s, value %d\n',kk, roiNames{kk}, roiCodes(kk)); 
        plot(avgROISig(kk,:),'DisplayName',sprintf('%s', roiNames{kk}));
    end
end
legend show; legend('FontSize',16);


%%
[u,s,d] = svd(avgROISig');


%% Save output
% reshapde the outputs to 4D 
avgGatedSignal = reshape(avgGatedSignal,[size(rsfmri,1:3),size(avgGatedSignal,2)]);
gatedSignal    = reshape(gatedSignal,   [size(rsfmri,1:3),size(gatedSignal,2)]);
gatedTime      = reshape(gatedTime,     [size(rsfmri,1:3),size(gatedTime,2)]);
% save avgGatedSignal into 4D NIFTI
ni_info = rsfmri_info;
ni_info.Datatype = 'double';
% avgGatedSignal
ni_info.ImageSize = size(avgGatedSignal);
niftiwrite(avgGatedSignal, fullfile(datadir,'avgGatedSignal'), ni_info);
% gatedSignal
ni_info.ImageSize = size(gatedSignal);
niftiwrite(gatedSignal, fullfile(datadir,'gatedSignal'), ni_info);
% gatedTime
ni_info.ImageSize = size(gatedTime);
niftiwrite(gatedTime, fullfile(datadir,'gatedTime'), ni_info);


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