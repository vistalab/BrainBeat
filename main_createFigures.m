

% This script generates figures from the manuscript titled:
%
%%%%% Measuring brain beats: cardiac-aligned fast fMRI signals %%%%%
% Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell
%

%% Check whether required paths are added

% check for SPM
try spm('Defaults','fmri')
catch ME
    switch ME.identifier
        case 'MATLAB:UndefinedFunction'
            warning('Add SPM12 to path (DO NOT USE genpath)');
    end
end

if ~exist('vistaRootPath.m','file')
    warning('Add vistasoft to Matlab path');
end

try bids.bids_matlab_version
catch ME
    switch ME.identifier
        case 'MATLAB:undefinedVarOrClass'
            warning('Add bids_matlab to Matlab path');
    end
end

[~,dDir] = bbPath;

%%

bb01_Method_Fig1

%%

bb02_OverlayR2Anat_Fig2

%%

bb02_OverlayR2MNI_Fig2J

%% 

bb03_CheckSignalSign_Fig3BC

%%

bb03_DKTatlasRender_Fig3D_Fig4

%%

bb05_SVD_AllSubs_Fig5_SuppFig4

%%

bb07_peaktimeMNI_Fig6

%%

bb08_peaktimeSubs_Fig7_Fig8

%%

bb09_restingStatefMRI_Fig9
