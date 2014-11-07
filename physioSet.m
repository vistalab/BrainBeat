function p = physioSet(p,param,val,varargin)
%Organize methods for setting physiological parameters.
%  
%   p = physioSet(p,param,val,varargin)
%
% Example:
%   p = physioCreate; p.ppg.srate
%   p = physioSet(p, 'ppg srate', 75); p.ppg.srate
%
% See also:  physioCreate, physioGet
%
%   
% DH/BW Vistasoft Team , 2014

if ~exist('val','var'), error('Value required'); end

param = mrvParamFormat(param);

switch param(1:3)
    case 'ppg'
        
        switch param(4:end)
            case 'name'
                p.ppg.name = val;
            case 'srate'
                p.ppg.srate = val;
            case 'rawdata'
                p.ppg.rawdata = val;
            case 'timeseries'
                p.ppg.timeseries = val;
            case 'scan_onset'
                p.ppg.scan_onset = val;
                
            otherwise
                error('Unknown PPG parameter %s\n',param);
        end
        
    case 'res'
        switch param(5:end)
            case 'name'
                p.resp.name = val;
            case 'srate'
                p.resp.srate = val;
            case 'rawdata'
                p.resp.rawdata = val; 
            case 'timeseries'
                p.resp.timeseries = val;
            case 'scan_onset'
                p.resp.scan_onset = val;
                
            otherwise
                error('Unknown RESP parameter %s\n',param);
        end
        
    otherwise
        error('Unknown physiological parameter %s\n',param);
        
end

end
