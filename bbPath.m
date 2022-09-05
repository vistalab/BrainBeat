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
addpath(rootPath);

addpath(fullfile(rootPath,'bbeat'));
addpath(genpath(fullfile(rootPath,'external')));
addpath(fullfile(rootPath,'physio'));
addpath(fullfile(rootPath,'roi_labeling'));
addpath(fullfile(rootPath,'scripts'));
addpath(fullfile(rootPath,'visualization'));


% make sure local data path is not added as code
rmpath(fullfile(rootPath,'local'));

% data directory
dDir = fullfile(rootPath,'local','BrainBeat');

return