
addpath(genpath('/share/wandell/users/dhermes/git/BrainBeat'))
addpath(genpath('/share/wandell/users/dhermes/git/vistasoft'))
addpath('/share/wandell/users/dhermes/git/spm12')


%% Base data directory
clear all 
close all
dDir = '/sni-storage/wandell/data/BrainBeat/data/';

%% Check the coregistration
s = 2;
s_info = bb_subs(s);
subj=s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

scan_nr = 1; % selects the functional scan, note that we're using the second scan here, with FA = 25|20, the other FA do not coregister well to the T1
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};

fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
ni = niftiRead(fmri);
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']),'acpcXform_new')

%% original data
% reshape into voxel X time
data_cat = reshape(ni.data,prod(ni.dim(1:3)),size(ni.data,4));
% get rid of first time points
data_cat = data_cat(:,5:end);

% calculate temporal snr
data_snr = mean(data_cat,2)./std(data_cat,[],2);
% reshape data_snr into 3D matrix again
data_snr = reshape(data_snr,ni.dim(1:3));

% calculate std
data_std = std(data_cat,[],2);
% reshape data_snr into 3D matrix again
data_std = reshape(data_std,ni.dim(1:3));

% calculate mean
data_mean = mean(data_cat,2);
% reshape data_snr into 3D matrix again
data_mean = reshape(data_mean,ni.dim(1:3));


%% for visualization:
curPos = [1,10,-30]; 
sliceThisDim = 2; 
% imDims=[-90 -120 -60; 90 130 90];
imDims=[-90 -120 -120; 90 130 90];

%%%% Overlay 1: snr
niFunc = ni;
niFunc.data = data_snr;
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform_new,sliceThisDim,imDims,curPos)
title('snr')

%%%% Overlay 2: std
niFunc = ni;
niFunc.data = data_std;
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform_new,sliceThisDim,imDims,curPos)
title('std')

%%%% Overlay 3: mean
niFunc = ni;
niFunc.data = data_mean;
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform_new,sliceThisDim,imDims,curPos)
title('mean')

%% now average across every 2 scans

% reshape into voxel X time
data_cat = reshape(ni.data,prod(ni.dim(1:3)),size(ni.data,4));
% get rid of first time points
data_cat = data_cat(:,6:end);

% create fake data to pretend data were sampled 2X slower, and every 2
% scans are averaged:
data_cat = (data_cat(:,1:2:end) + data_cat(:,2:2:end))./2;

% calculate temporal snr
data_snr = mean(data_cat,2)./std(data_cat,[],2);
% reshape data_snr into 3D matrix again
data_snr = reshape(data_snr,ni.dim(1:3));

% calculate std
data_std = std(data_cat,[],2);
% reshape data_snr into 3D matrix again
data_std = reshape(data_std,ni.dim(1:3));

% calculate mean
data_mean = mean(data_cat,2);
% reshape data_snr into 3D matrix again
data_mean = reshape(data_mean,ni.dim(1:3));


%% for visualization:
curPos = [1,10,-30]; 
sliceThisDim = 2; 
% imDims=[-90 -120 -60; 90 130 90];
imDims=[-90 -120 -120; 90 130 90];

%%%% Overlay 1: snr
niFunc = ni;
niFunc.data = data_snr;
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform_new,sliceThisDim,imDims,curPos)
title('snr')

%%%% Overlay 2: std
niFunc = ni;
niFunc.data = data_std;
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform_new,sliceThisDim,imDims,curPos)
title('std')

%%%% Overlay 3: mean
niFunc = ni;
niFunc.data = data_mean;
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform_new,sliceThisDim,imDims,curPos)
title('mean')

