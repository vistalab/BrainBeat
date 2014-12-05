function p = physioCreate(varargin)
% Create a physiological data structure for heart and respiration
%
%  p = physioCreate(varargin)
%
% The heart beat is estimated with a photoplethysmogram
% The respiration is estimated with a breathing belt.
%
% Example:
%    p = physioCreate;
%    p = physioCreate('filename',fname);
%     
%    p = physioCreate('ppg srate',50);
%
%    p = physioCreate('resp srate',33)
%
% See also:  bbGet
%
% DH/BW Vistasoft team, 2014

ppg = ppgInit;
p.ppg = ppg;

resp = respInit;
p.resp = resp;

p.filename = 'filename';

if ~isempty(varargin)
    if ~mod(length(varargin),2)  % Must be param/val pairs
        if strcmp(varargin{1},'filename')
            p.filename = varargin{2};
            % Do that bbGet thing in here to read the file data
            % Figure that the file is right.
            return
        else
            % Overwrite default arguments with the input arguments
            for ii=1:2:length(varargin)
                p = physioSet(p,varargin{ii},varargin{ii+1});
            end
        end
    end
end

end


%%%%
function ppg = ppgInit
% Physiological photoplethethysmogram data structure
%
%
% Time scale is always at the CNI sampling rate
%
%
ppg.name = 'PPG';
ppg.srate = 100;        % Number of samples per second

ppg.rawdata    = [];    % Values in the file
ppg.timeseries = [];
ppg.scan_onset = [];    % Times when the scan acq begins

end

%%%%
function resp = respInit

resp.name = 'RESP';
resp.srate = 25;

resp.rawdata    = [];    % Values in the file
resp.timeseries = [];
resp.scan_onset = [];    % Times when the scan acq begins

end
