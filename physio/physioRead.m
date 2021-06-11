function [rawdata, syncdata] = physioRead(fname, dt, scanDuration, dataType, scanStart)
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

rawdata = readmatrix(fname,'leadingDelimitersRule','ignore');

if strcmp(dataType, 'wave')
    syncdata = rawdata(end-round(scanDuration/dt)+1 : end);
elseif strcmp(dataType, 'trig')
    if ~exist('scanStart', 'var')
        error('Must provide scan start time to align the triggers \n');
    end
    syncdata = rawdata(rawdata * dt > scanStart) * dt - scanStart; % trigger events in milliseconds
end
    
end