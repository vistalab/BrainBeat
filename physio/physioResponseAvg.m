function [avgGatedSignal, interpWindows] = physioResponseAvg(gatedSignal, gatedTime, varargin)
% Averaging the cardiac gated signal to a specific time resolution
%
% Synopsis
%   [avgGatedSignal, interpTimes] = physioResponseAvg(gatedSignal, gatedTime, varargin)
%
% Inputs
%    gatedSignal    - original cardiac gated signal 
%    gatedTimes     - original cardiac delay time
%
% Optional key/val pairs
%    TR             - in ms, default = 2000
%    timeWindow     - in ms, how we temporally average the gated signal, default = 50
%    timeResolution - in ms, time resolution of the interpolated gated times, default = 10
%
% Outputs
%    avgGatedSignal - T2* values gated with respect to the heart beats
%                     Size should be (size(gatedSignal,1), nTimeSamples)
%    interpWindows  - time grid over which the gated signal is averaged
%
% See also
%   physioResponse
%

%% Parameters

p = inputParser;
p.addRequired('gatedSignal',@ismatrix);
p.addRequired('gatedTime',@isnumeric);
p.addParameter('TR',2000,@(x)(isnumeric(x) && (x == round(x))));
p.addParameter('timeWindow',50,@(x)(isnumeric(x) && (x == round(x))));
p.addParameter('timeResolution',10,@(x)(isnumeric(x) && (x == round(x))));

parse(p, gatedSignal, gatedTime, varargin{:});

TR           = p.Results.TR;
timeWindow   = p.Results.timeWindow;
deltaT       = p.Results.timeResolution;  % ms, resoluion of the time grid to interpolate

% The user may want the signal without the averaging
interpTimes = 0:deltaT:(TR-timeWindow);  % deltaT is the time steps
interpWindows = interpTimes + timeWindow/2;  % center of the time grids

% reshape gatedSignal into a 2D matrix in order to iterate over all voxels 
nVoxels = size(gatedSignal,1);
avgGatedSignal = zeros(nVoxels,length(interpTimes));

% For each voxel in the volume
for vv = 1:nVoxels
   
    % Calculate mean of the gated signal over the time windows
    for tt = 1:length(interpTimes)
        signalInWindow = gatedSignal(vv,(interpTimes(tt) <= gatedTime(vv,:) & gatedTime(vv,:) < interpTimes(tt) + timeWindow));
        avgGatedSignal(vv,tt) = mean(signalInWindow);
    end
    
end

end



