function [gatedSignal, gatedTimes] = physioResponse(hbeats, tseries, varargin)
% Gate rs time series to heart beats to look for response
%
% Synopsis
%   [gatedSignal, gatedTimes] = physioResponse(hbeats, tseries, varargin)
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
%
% Outputs
%    gatedSignal    - original cardiac gated signal 
%    gatedTimes     - original cardiac delay time
%
% See also
%   physioGet/Set/Create/Detrend
%

%% Parameters

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

parse(p, hbeats, tseries, varargin{:});

TR           = p.Results.TR;
slice_timing = p.Results.slice_timing;
roi          = p.Results.roi;
dummy_cycles = p.Results.dummy_cycles;

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

% Only keep the heart beats after the dummy cycles
hbeats = hbeats(hbeats > TR * dummy_cycles);
nBeats   = numel(hbeats);

% Create a vector of (gatedTime, T2*) for all voxels
gatedTimes  = zeros(nPoints,nBeats);
gatedSignal = zeros(nPoints,nBeats);

% Relate the row/col to the NIfTI format and make sure row/col is right
% Confirmed x is row, y is col in NIfTI. Matlab imagesc plots the figure with X/Y
% axes swapped.
row = roi(:,1); col = roi(:,2); sl = roi(:,3);

%% Find the fmri volume closest to each ppg trigger

% For each roi point
for rr = 1:nPoints

    x = row(rr); y = col(rr); z = sl(rr);
    % For this voxel, here are the slice times
    tt = ttSlices(:,z);
    
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
            gatedTimes(rr,trig) = tseries_vol_t - hbeats(trig);
            % This is the T2* signal at that time.
            gatedSignal(rr,trig) = tseries(x, y , z, idx);
        else
            gatedTimes(rr,trig) = NaN;
            gatedSignal(rr,trig) = NaN;
        end
        
    end
        
end

end



