function [avgGatedSignal,gatedSignal] = physioResponse(hbeats, tseries, TR, slice_timing, roi, dummy_cycles,avgT)
% Gate rs time series to heart beats to look for response
%
% Synopsis
%   [gatedSignal, avgGatedSignal] = physioResponse(hbeats, tseries, TR, nvols, slice_timing, roi, dummy_cycles)
%
% Inputs
%    hbeats  -  Vector of heart beat times in milliseconds
%    tseries -  4D T2* data set from a NIfTI
%    params  -  TR in ms,
%               slice_timing in ms,
%               roi - spatial coordinates of an ROI in (x,y,sl) format
%                     N x 3 matrix
%               dummy_cycles - Skip these nvols
%               timeWindow - milliseconds, how we temporally average the
%                  gated signal
%
% Optional key/val pairs
%
% Outputs
%   gatedSignal -  T2* values gated with respect to the heart beats
%                  Size should be (nTimeSamples x nROIPoints)
%
% See also
%   physioGet/Set/Create
%

%% Parameters

% physioResponse(hbeats,tseries,'dummy_cycle',3,'avgT',8);

% p = inputParser;
% p.addRequired('hbeats',@ismatrix);
% p.addRequired('tseries',@(x)(ismatrix(x) && (ndims(x) == 4)));
%
% p.addParameter('dummy_cycles',3,@isinteger);
% ...
%
% p.parse(

if ~exist('dummy_cycles', 'var')
    dummy_cycles = 3;
end
avgT = 8; 
deltaT  = 10;  % ms

nPoints = size(roi,1);
avgT    = 8;   % We average over this many deltaT for the avgGatedSignal

% nRows    = size(tseries,1);
% nCols    = size(tseries,2);
nSlices  = size(tseries,3);
nVolumes = size(tseries,4);

% The nominal timing of the volumes, not correcting for slice timing
t_volumes = (0 : (nVolumes - 1)) * TR;
ttSlices = zeros(length(dummy_cycles:nVolumes),nSlices);
for ss = 1:nSlices
    ttSlices(:,ss) = t_volumes(dummy_cycles:end) + slice_timing(ss);
end

% We return the average of the gated signal by default
avgGatedSignal = zeros(nvols,nPoints);

% The user may want the signal without the averaging
interpTimes = deltaT:deltaT:TR;  % deltaT is the time tsteps

% Only keep the heart beats after the dummy cycles
hbeats = hbeats(hbeats > TR * (dummy_cycles-1));
nBeats   = numel(hbeats);

% There will be 
gatedSignal = zeros(nBeats,nPoints);

%% Remove a 2nd order polynomial from the time series

tseries_reg = physioDetrend(tseries, dummy_cycles);

%% Find the fmri volume closest to each ppg trigger

% For each roi point
for rr = 1:nPoints
    % Relate the row/col to the NIfTI format and make sure row/col is right
    col = roi(rr,1); row = roi(rr,2); sl = roi(rr,3);

    % For this voxel, here are the slice times
    tt = ttSlices(:,sl);
    
    % Create a vector of (gatedTime, T2*) for this voxel
    gatedTime   = zeros(nBeats);
    sig_hbeats = zeros(nBeats);
    for trig = 1:nBeats
        % For each heart beat, find the first slice time after that beat.
        % In principle, we could take all the slices within some temporal
        % window of the heart beat.  But for a 2sec TR, and heart beats
        % about 1 sec, we decided to just take the nearest volume.
        idx = find(tt > hbeats(trig), 1, 'first');
        
        % This is the time in milliseconds of that slice
        tseries_vol_t    = tt(idx);
        
        % This is the time since the heart beat
        gatedTime(trig)    = tseries_vol_t - hbeats(trig);
        
        % This is the T2* signal at that time.
        sig_hbeats(trig) = tseries_reg(col, row , sl, idx);
        
    end
        
    % Sort the data by gated time
    [sortedGatedTimes, idx] = sort(gatedTime);
    gatedSignal(:,rr) = sig_hbeats(idx);
    
    % Smooth and interpolate the time base for the returned signal
    
    % Make sure there are no precise overlaps to stop interp1 from
    % complaining
    uniqueTimes = sortedGatedTimes + rand(size(sortedGatedTimes));
    sig_interp = interp1(uniqueTimes, gatedSignal(:,rr), interpTimes);
    
    % Just a box average over some number of periods.
    % avgT * deltaT is the averaging time
    ker = ones(1, avgT);
    avgGatedSignal(:,rr) = conv(sig_interp, ker, 'same');
    
end

end



