function val = physioGet(phy,param,varargin)
% Get for physio parameters
%
%   val = physioGet(phy,param,varargin)
%
% Input:
%   phy:   Physiological data structure, including ppg and resp
%   param: Parameter to retrieve
% 
% Output:
%   val:  return value
%
% Examples:
%   We should create a test physio file that we can load easily and then
%   write unit tests for the gets in this function.
%
% See also: physioSet, physioCreate
% 
% DH/Vistasoft Team Copyright Stanford 2014

if ~exist('phy','var') || isempty(phy), error('phy required'); end
if ~exist('param','var') || isempty(param), error('param required'); end

param = mrvParamFormat(param);

% If there is a filename, return it.
if strcmp(param,'filename')
    if isfield(phy,'filename'), val = phy.filename; return; 
    else error('No filename field');
    end
end

% Which type of data are we getting?
respFlag = false; ppgFlag = false;
if     strcmp(param(1:4),'resp'), respFlag = true; param=param(5:end); 
elseif strcmp(param(1:3),'ppg'),  ppgFlag = true;  param = param(4:end);
end

% Respiratory data
if respFlag
    
    switch param
        case 'data'
            val = phy.resp.data;
        case 'peaks'
            % physioGet(phy,'resp peaks',[interval])
            if isempty(varargin), interval = 2;
            else interval = varargin{1};
            end
            signal = physioGet(phy,'resp data');
            srate  = physioGet(phy,'resp srate');

            % set minimum inter-heartbeat interval in seconds
            val = physioPeaks(signal,interval,srate);
            val = val / srate;

        case 'srate'
            % Number of samples per second
            val = phy.resp.srate;
            
        case 'sampletimes'
            % Could adjust for sec or ms or other units.  Ask BW.
            % seconds/persample
            val = (1:physioGet(phy,'resp n samples'))/physioGet(phy,'resp srate');
            
        case 'nsamples'
            val = length(phy.resp.data);
            
        otherwise
            error('Unknown parameter %s\n',param);
    end
    
% Pulse data from the plethysmograph
elseif ppgFlag
    
    switch param
        case 'data'
            % Raw data
            val = phy.ppg.data;
        case 'peaks'
            % Returns the time of the peaks in secs
            if isempty(varargin), interval = 0.7;
            else interval = varargin{1};
            end
            signal = physioGet(phy,'ppg data');
            srate  = physioGet(phy,'ppg srate');

            
            % set minimum inter-heartbeat interval in seconds
            val = physioPeaks(signal,interval,srate);
            val = val / srate;
        case 'srate'
            % Samples per second
            val = phy.ppg.srate;
            
        case 'sampletimes'
            % Could adjust for sec or ms or other units.  Ask BW.
            val = (1:physioGet(phy,'ppg n samples'))/physioGet(phy,'ppg srate');
            
        case 'nsamples'
            val = length(phy.ppg.data);
            
            
        otherwise
            error('Unknown parameter %s\n',param);
    end
    
end

end

% ----
function onsets = physioPeaks(signal,interval,srate)
% Calculate the peaks for either the PPG or RESP data

% low-pass filter 
% not necessary, keep in for lower quality data?
band = 5;
Rp   = 3; Rs=60; % third order Butterworth
high_p =  band(1)*2/srate;
high_s = (band(1)+20)*2/srate;

[n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
[bf_b,bf_a] = butter(n_band,wn_band,'low');
band_sig    = filtfilt(bf_b,bf_a,signal);

% detect peak sample positions
[~,onsets] = findpeaks(band_sig,'minpeakdistance',interval*srate);

end




        
