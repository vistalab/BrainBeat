function [rootPath,dDir] = bbPath()
% return root of the BrainBeat code and data directories
%
% [rootPath dDir] = bbPath;
%
% This function MUST reside in the directory at the base of the BrainBeat
% directory structure
%
% Wandell Copyright Vistasoft Team, 2013

% code directory
rootPath = fileparts(which(mfilename));

% add code path
addpath(genpath(rootPath));

% make sure local data path is not added as code
rmpath(fullfile(rootPath,'local'));

% data directory
dDir = fullfile(rootPath,'local','BrainBeat');

return