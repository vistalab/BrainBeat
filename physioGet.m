function val = physioGet(phy,param,varargin)
% Get for physio parameters
%
%   val = physioGet(phy,param,varargin)
%
% Examples:
%
% See also: physioSet, physioCreate
% 
% DH/Vistasoft Team Copyright Stanford 2014

if ~exist('phy','var') || isempty(phy), error('phy required'); end
if ~exist('param','var') || isempty(param), error('param required'); end

param = mrvParamFormat(param);

% We should switch on ppg and resp parameters in some sensible way here.
% It could be that the names are always
%   physioGet(p,'resp/ppg param') and we detect which at the beginning
%

switch param
    case 'filename'
        val = phy.filename;
    otherwise
end

end
