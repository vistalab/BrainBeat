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
%   val:  Return value
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
    obj = phy.resp;
    dataType='resp';
elseif strcmp(param(1:3),'ppg'),  ppgFlag = true;  param = param(4:end);
    obj = phy.ppg;
    dataType='ppg';
end


switch param
    case 'data'
        val = obj.data;
    case 'peaks'
        % physioGet(phy,'resp peaks',[interval])
        
        signal = physioGet(phy,[dataType 'data']);
        srate  = physioGet(phy,[dataType 'srate']);
        
        if isempty(varargin)
        % we could use preset minimum interval
%             if respFlag
%                 interval = 2; 
%             elseif ppgFlag 
%                 interval = .7; 
%             end
        % but let's get it from the data:
            % find the peaks in the autocorrelation function:
            [peaks_ac,peaks_ac_i] = findpeaks(acf(obj.data,500));
            % the first maximum peak (zero is not included) is the first
            % autocorrelation time
            [~,max_peaks_ac_i]=max(peaks_ac);
            % set the interval at 70% of the ppg/resp rate:
            interval = .7 * peaks_ac_i(max_peaks_ac_i)/srate;
        else interval = varargin{1};
        end
        
        % set minimum inter-heartbeat interval in seconds
        val = physioPeaks(signal,interval,srate);
        val = val / srate;

    case 'srate'
        % Number of samples per second
        val = obj.srate;

    case 'sampletimes'
        % Could adjust for sec or ms or other units.  Ask BW.
        % seconds/persample
        val = (1:physioGet(phy,[dataType ' n samples']))/physioGet(phy,[dataType ' srate']);

    case 'nsamples'
        val = length(obj.data);

    case 'rate' % peaks per second
        objPeaks = physioGet(phy,[dataType ' peaks']);
        totalTime = length(obj.data)./physioGet(phy,[dataType ' srate']); % in secs
        val = length(objPeaks) / totalTime; 
            
    case 'ppgcurve' % one response curve, starting at peak -0.5 sec
        objPeaks = physioGet(phy,[dataType ' peaks']);
        objRate = physioGet(phy,[dataType ' rate']);
        srate = physioGet(phy,[dataType ' srate']);
        signal = physioGet(phy,[dataType 'data']);
        
%         epoch_length = ceil(srate/objRate); % this would be one response curve
        epoch_length = floor(srate*2.5); % 2.5 sec 
        objPeaks(objPeaks*srate-floor(srate*.5)<1)=[]; % get rid of first Peak, may nog capture full response
        objPeaks(objPeaks*srate+epoch_length>length(signal))=[]; % get rid of last Peak, may nog capture full response
        phys_epochs1 = zeros(length(objPeaks),epoch_length);
        for k = 1:length(objPeaks)
            sample_ppg=round(objPeaks(k)*srate);
            phys_epochs1(k,:) = signal(sample_ppg-floor(srate*.5):sample_ppg+epoch_length-floor(srate*.5)-1);
        end
        
        % get average ppg response
        val=mean(phys_epochs1,1);
        
    case 'ppgtcurve' % time for one response curve
        
        objRate = physioGet(phy,[dataType ' rate']);
        srate = physioGet(phy,[dataType ' srate']);       
        epoch_length = floor(srate*2.5);
        tcurve = [1:epoch_length]/srate -.5;
        
        val = tcurve;
    
    case 'respcurve' % one response curve, starting at peak -0.5 sec
        objPeaks = physioGet(phy,[dataType ' peaks']);
        objRate = physioGet(phy,[dataType ' rate']);
        srate = physioGet(phy,[dataType ' srate']);
        signal = physioGet(phy,[dataType 'data']);
        
        epoch_length = ceil(srate/objRate); % this would be one response curve
        objPeaks(objPeaks*srate+epoch_length>length(signal))=[]; % get rid of last Peak, may nog capture full response
        phys_epochs1 = zeros(length(objPeaks),epoch_length);
        for k = 1:length(objPeaks)
            sample_ppg=round(objPeaks(k)*srate);
            phys_epochs1(k,:) = signal(sample_ppg:sample_ppg+epoch_length-1);
        end
        
        % get average ppg response
        val=mean(phys_epochs1,1);
        
    case 'resptcurve' % time for one response curve
        
        objRate = physioGet(phy,[dataType ' rate']);
        srate = physioGet(phy,[dataType ' srate']);       
        epoch_length = ceil(srate/objRate); % .5 sec + one response curve
        tcurve = [1:epoch_length]/srate;
        
        val = tcurve;
        
    otherwise
        error('Unknown parameter %s\n',param);
end


end

% ----
function onsets = physioPeaks(signal,interval,srate)
% Calculate the peaks for either the PPG or RESP data

% low-pass filter 
% not necessary, keep in for lower quality data?
band = 5;
Rp   = 3; Rs = 60; % third order Butterworth
high_p =  band(1)*2/srate;
delta = 0.001*2/srate;
high_s = min(1-delta,high_p+0.1);

[n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
[bf_b,bf_a] = butter(n_band,wn_band,'low');
band_sig    = filtfilt(bf_b,bf_a,signal);

% detect peak sample positions
[~,onsets] = findpeaks(band_sig,'minpeakdistance',round(interval*srate));

end




        
