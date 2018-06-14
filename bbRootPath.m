function rootPath=bbPath()
% return root of the BrainBeat directory
%
% rootPath = bbPath;
%
% This function MUST reside in the directory at the base of the BrainBeat
% directory structure
%
% Wandell Copyright Vistasoft Team, 2013

rootPath = fileparts(which(mfilename));

return