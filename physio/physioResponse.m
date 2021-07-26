function [avgGatedSignal, interpWindows, gatedSignal, gatedTimes] = physioResponse(hbeats, tseries, varargin)
% Gate rs time series to heart beats to look for response
%
% Synopsis
%   [avgGatedSignal, interpTimes, sortedGatedSignal, sortedGatedTimes] = physioResponse(hbeats, tseries, varargin)
%
% Inputs
%    hbeats         - Vector of heart beat times in milliseconds
%    tseries        - 4D T2* data set after detrending for a 2nd order polynomial

% Optional key/val pairs
%    TR             - in ms, default = 2000
%    slice_timing   - in ms, default is calculated with 2000ms TR
%    roi            - spatial coordinates of an ROI in (x,y,sl) format 
%                      N x 3 matrix, default is the center of volume
%    dummy_cycles   - Skip these nvols, default = 3
%    timeWindow     - in ms, how we temporally average the gated signal, default = 50
%    timeResolution - in ms, time resolution of the interpolated gated times, default = 10
%
% Outputs
%    avgGatedSignal - T2* values gated with respect to the heart beats
%                     Size should be (nTimeSamples x nROIPoints)
%    interpWindows  - time grid over which the gated signal is averaged
%    gatedSignal    - original cardiac gated signal 
%    gatedTimes     - original cardiac delay time
%
% See also
%   physioGet/Set/Create/Detrend
%

%% Parameters

% physioResponse(hbeats,tseries,'dummy_cycle',3,'avgT',8);

defaultTR = 2000;
[~,defaultSliceOrder] = sort([1:2:size(tseries,3), 2:2:size(tseries,3)]);  % assuming non-SMS interleaved acquisition
defaultSliceTiming =  defaultSliceOrder/ size(tseries,3) * defaultTR;
defaultROI = floor(size(tseries, 1:3)/2);       % center of the volume

p = inputParser;
p.addRequired('hbeats',@ismatrix);
p.addRequired('tseries',@(x)(isnumeric(x) && (ndims(x) == 4)));
p.addParameter('TR',2000,@(x)(isnumeric(x) && (x == round(x))));
p.addParameter('slice_timing',defaultSliceTiming,@isnumeric);
p.addParameter('roi',defaultROI,@ismatrix);
p.addParameter('dummy_cycles',3,@(x)(isnumeric(x) && (x == round(x))));
p.addParameter('timeWindow',50,@(x)(isnumeric(x) && (x == round(x))));
p.addParameter('timeResolution',10,@(x)(isnumeric(x) && (x == round(x))));

parse(p, hbeats, tseries, varargin{:});

TR           = p.Results.TR;
slice_timing = p.Results.slice_timing;
roi          = p.Results.roi;
dummy_cycles = p.Results.dummy_cycles;
timeWindow   = p.Results.timeWindow;
deltaT       = p.Results.timeResolution;  % ms, resoluion of the time grid to interpolate
avgT         = timeWindow/deltaT;   % We average over this many deltaT for the avgGatedSignal

nPoints = size(roi,1);

% nRows    = size(tseries,1);
% nCols    = size(tseries,2);
nSlices  = size(tseries,3);
nVolumes = size(tseries,4);

% The nominal timing of the volumes, not correcting for slice timing
t_volumes = ((0:(nVolumes-1)) + dummy_cycles) * TR;
ttSlices = zeros(nVolumes,nSlices);
for ss = 1:nSlices
    ttSlices(:,ss) = t_volumes + slice_timing(ss);
end

% The user may want the signal without the averaging
interpTimes = 0:deltaT:(TR-timeWindow);  % deltaT is the time steps
interpWindows = interpTimes + timeWindow/2;  % center of the time grids

% We return the average of the gated signal by default
avgGatedSignal = zeros(length(interpTimes),nPoints);

% Only keep the heart beats after the dummy cycles
hbeats = hbeats(hbeats > TR * dummy_cycles);
nBeats   = numel(hbeats);

% There will be nBeats T2* values and delay times for each roi point
% sortedGatedSignal = zeros(nBeats,nPoints);
% sortedGatedTimes  = zeros(nBeats,nPoints);

%% Find the fmri volume closest to each ppg trigger

% For each roi point
for rr = 1:nPoints
    % Relate the row/col to the NIfTI format and make sure row/col is right
    % Confirmed x is row, y is col in NIfTI. Matlab imagesc plots the figure with X/Y
    % axes swapped.
    row = roi(rr,1); col = roi(rr,2); sl = roi(rr,3);

    % For this voxel, here are the slice times
    tt = ttSlices(:,sl);
    
    % Create a vector of (gatedTime, T2*) for this voxel
    gatedTimes  = zeros(1, nBeats);
    gatedSignal = zeros(1, nBeats);
    for trig = 1:nBeats
        % For each heart beat, find the first slice time after that beat.
        % In principle, we could take all the slices within some temporal
        % window of the heart beat.  But for a 2sec TR, and heart beats
        % about 1 sec, we decided to just take the nearest volume.
        idx = find(tt > hbeats(trig), 1, 'first');

        if ~isempty(idx)
            % This is the time in milliseconds of that slice
            tseries_vol_t    = tt(idx);
            % This is the time since the heart beat
            gatedTimes(trig) = tseries_vol_t - hbeats(trig);
            % This is the T2* signal at that time.
            gatedSignal(trig) = tseries(row, col , sl, idx);
        else
            gatedTimes(trig) = NaN;
            gatedSignal(trig) = NaN;
        end
        
    end
        
    % Sort the data by gated time
    % [sortedGatedTimes(:,rr), idx] = sort(gatedTimes);
    % sortedGatedSignal(:,rr) = gatedSignal(idx);
    
    % Smooth and interpolate the time base for the returned signal
    
    % Make sure there are no precise overlaps to stop interp1 from
    % complaining
    % uniqueTimes = sortedGatedTimes + rand(size(sortedGatedTimes))*0.001;
    % sig_interp = interp1(uniqueTimes, sortedGatedSignal(:,rr), interpTimes);
    % figure;plot(interpTimes,sig_interp,'x')
    
    % Just a box average over some number of periods.
    % avgT * deltaT is the averaging time
    % ker = ones(1, avgT);
    % avgGatedSignal(:,rr) = conv(sig_interp, ker, 'same') / avgT;
    
    % Calculate mean of the gated signal over the time windows
    for tt = 1:length(interpTimes)
        signalInWindow = gatedSignal((interpTimes(tt) <= gatedTimes & gatedTimes < interpTimes(tt) + timeWindow));
        avgGatedSignal(tt, rr) = mean(signalInWindow);
    end
    
end

end



