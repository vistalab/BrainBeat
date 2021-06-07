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

%% Here are the scitran commands for finding and downloading the data

% see bbResponse2physio in BrainBeats

% Let's try with scitran for a while.
st = scitran('cni');

% Or use Flywheel MATLAB SDK. fw = flywheel.Flywheel(apikey)
% addpath('~/Documents/MATLAB/'); fw = setup_flywheel;

%% Download files from Flywheel using scitran

% General information so we can find the physio data
group = 'leanew1';
project_label = 'iris';
subj_label = 'IRIS021';
session_label = '22776';
acq_label = 'resting_2';

% Several ways to find data

% Using lookup 
proj = st.lookup(fullfile(group,project_label));

% Using the FW methods and stSelect
sessions = proj.sessions();
thisSession = stSelect(sessions,'label',session_label);

% Using find
% thisSession = proj.sessions.find(sprintf('label="%s"',session_label))

% Examples of lookup methods that include the subject
lookup = fullfile(group,project_label,subj_label,session_label,acq_label);
thisAcq = st.lookup(lookup);

% See the file names and download and unzip
stPrint(thisAcq.files,'name')

chdir(fullfile(bbPath,'local'));
gephysioZip = stSelect(thisAcq.files,'name','gephysio','contains',true);
gephysioFile = fullfile(bbPath,'local','geophysio.zip');
gephysioZip{1}.download(gephysioFile);
unzip(gephysioFile);

%{
% Get the physio data and acquisition information from Flywheel if SDK exists 
if exist('fw', 'var')
    group = 'cni';
    project_label = 'iris';
    session_date = '2021-02-20';
    acquisition = 'resting_2';
    % acq = st.lookup('leanew1/iris/2021-02-20/resting_2')
    proj = fw.projects.findFirst(['label=' project_label]);
    sess = proj.sessions.find(['created>' session_date]);
    sess = sess{1};
    acq = sess.acquisitions.find(['label=' acquisition]);
    acq = acq{1};
    files = acq.files;
    for f=1:numel(files)
        if strcmp(files{f}.type,'gephysio')  % download the *.gephysio.zip and unzip to gephysio/data
            acq.downloadFile(files{f}.name, fullfile(dataRootPath, files{f}.name));
            unzip(fullfile(dataRootPath, files{f}.name), dataRootPath);
            [~,subfolder,~] = fileparts(files{f}.name);
            datadir = fullfile(dataRootPath, subfolder);
        end
        if strcmp(files{f}.type,'dicom')   % get one dicom image from the dicom archive
            zipInfo = acq.getFileZipInfo(files{f}.name);
            dcmCount = numel(zipInfo.members);
            dcmName = zipInfo.members{1}.path;
            acq.downloadFileZipMember(files{f}.name, dcmName, fullfile(datadir, '1.dcm'));
            dcm = dicominfo(fullfile(datadir, '1.dcm'));
            % Get number of volumes, TR from the dicom header. Assuming this is a complete scan
            nvols = dcm.NumberOfTemporalPositions;  
            TR = dcm.RepetitionTime;
            scan_duration = TR * nvols;
       end
    end
else   % or manually prepare data
    datadir = fullfile(dataRootPath, 'gephysio1');
    TR = 2000;
    nvols = 240;
    scan_duration = TR * nvols;
end
%}

%% pre-set some parameters
ppg_sample = 10;        % PPG data samples at 10ms
resp_sample = 40;       % Respiratory data samples  at 40 ms

plot_window = 10000;    % Duration of physio data to plot for inspection (10 seconds)


%% Pulsimetry data
f = dir(fullfile(datadir,'PPGData*'));
foo = readmatrix(fullfile(f.folder, f.name));
ppg = foo(:,2);
f = dir(fullfile(datadir,'PPGTrig*'));
ppg_trig = readmatrix(fullfile(f.folder, f.name)); % These are times that the PPG detects as an event

% Convert to millisecond. 
ppg_t = ((1:numel(ppg))-1)*ppg_sample;
ppg_trig = ppg_trig*ppg_sample;

% Choose a range of times to plot
start = ppg_t(end) - plot_window;
stop  = ppg_t(end);

ppg_r = logical(start < ppg_t) & logical(ppg_t < stop);
ppg_mx = max(ppg(ppg_r)); ppg_mn = min(ppg(ppg_r));
ppg_lst = logical(start < ppg_trig) & logical (ppg_trig < stop);

figure; subplot(2,1,1);
plot(ppg_t(ppg_r),ppg(ppg_r)); hold on;
plot(ppg_trig(ppg_lst),ones(size(ppg_trig(ppg_lst)))*ppg_mx,'ko');
xlabel('Time (ms)'); ylabel('PPG recorded data'); set(gca,'ylim',[ppg_mn ppg_mx], 'xlim',[start, stop]);

%%  Respiration data
f = dir(fullfile(datadir,'RESPData*'));
foo = readmatrix(fullfile(f.folder, f.name));
resp = foo(:,2);
f = dir(fullfile(datadir,'RESPTrig*'));
resp_trig = readmatrix(fullfile(f.folder, f.name));

resp_t = ((1:numel(resp))-1)*resp_sample;
resp_trig = resp_trig*resp_sample;
resp_r = logical(start < resp_t) & logical(resp_t < stop);
resp_lst = logical(start < resp_trig) & logical (resp_trig < stop);
resp_mx = max(resp(resp_r)); resp_mn = min(resp(resp_r));

subplot(2,1,2);
plot(resp_t(resp_r),resp(resp_r)); hold on;
plot(resp_trig(resp_lst),ones(size(resp_trig(resp_lst)))*resp_mx,'ko');
xlabel('Time (ms)'); ylabel('RESP recorded data'); set(gca,'ylim',[resp_mn resp_mx],'xlim',[start, stop]);

%% Write out simple summary of the data
s.total_recording_time = ppg_t(end);    % in milliseconds
s.scan_duration = scan_duration;        % in milliseconds
s.scan_start = s.total_recording_time - scan_duration;
s.ppg_data_during_scan  = ppg((ppg_t(end) - scan_duration)/ppg_sample : end);     % one data point per 10ms
s.ppg_trig_during_scan  = ppg_trig(ppg_trig > s.scan_start) - s.scan_start;       % time position of events in milliseconds
s.resp_data_during_scan = resp((resp_t(end) - scan_duration)/resp_sample : end);  % one data point per 40ms
s.resp_trig_during_scan = resp_trig(resp_trig > s.scan_start) - s.scan_start;     % time position of events in milliseconds

fid = fopen(fullfile(datadir, 'physio.json'), 'w');
fprintf(fid, jsonencode(s, 'PrettyPrint', true));
fclose(fid);

%% END