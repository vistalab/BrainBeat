function [rawdata, syncdata, t_raw, t_sync, scanStart] = physioRead(fname, dt, scanDuration, dataType, scanStart)
% Read the physio files and figure out the synchronization
%
% Input:
%   fname           file name of the physio data
%   dt              time resolution of the physio data
%   scanDuration    duration of the MRI data acquisition
%   dataType        type of the GE physio data, waveform (*Data) or trigger (*Trig)
%   scanStart       starting time of the scan, relative to the beginning of
%                   the physio data recording. Optional for waveform data, 
%                   required for trigger data
%
% Output:
%   rawdata         orignal data in the physio file
%   syncdata        data after synchronizing with MRI acquisition
%   t_raw           time samples in the raw physio waveform. Ignored for trigger data
%   t_sync          time samples in the synced physio waveform. Ignored for trigger data
%   scanStart       starting time of the scan, calculated from the waveform

rawdata = readmatrix(fname,'leadingDelimitersRule','ignore');

if strcmp(dataType, 'wave')
    syncdata = rawdata(end-round(scanDuration/dt)+1 : end);
    scanStart = length(rawdata)*dt - scanDuration; % in milliseconds
    t_sync = linspace(dt, scanDuration, numel(syncdata));
    t_raw  = linspace(dt, length(rawdata)*dt, length(rawdata));
elseif strcmp(dataType, 'trig')
    if ~exist('scanStart', 'var')
        error('Must provide scan start time to align the triggers \n');
    end
    syncdata = rawdata(rawdata * dt > scanStart) * dt - scanStart; % trigger events in milliseconds
    
end

if ~exist('t_sync', 'var'), t_sync = []; end
if ~exist('t_raw', 'var'),  t_raw  = []; end
  
end