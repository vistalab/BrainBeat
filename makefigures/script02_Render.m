clear all
close all

%% Make a rendering 

% mris_convert lh.pial lh.pial.gii

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

s_nr = 2;
s_info = bb_subs(s_nr);
subj=s_info.subj;

%% Get the Flywheel obj file and make it a gifti:

mri_orig = fullfile(dDir,subj,'freesurfer','mri','orig.mgz');
fw_obj(1).name = fullfile(dDir,subj,'freesurfer','FlyWheelAnalyses','lh.pial.obj');
fw_obj(2).name = fullfile(dDir,subj,'freesurfer','FlyWheelAnalyses','rh.pial.obj');
fw_obj_hemi = {'l','r'};

% Get transformation matrix from freesurfer obj to original MRI:
orig = MRIread(mri_orig);
Torig = orig.tkrvox2ras;
Norig = orig.vox2ras;
freeSurfer2T1 = Norig*inv(Torig);

for kk = 1:length(fw_obj)
    [vertex,face] = read_obj(fw_obj(kk).name);
    g.vertices = vertex';
    g.faces = face';
    g.mat = eye(4,4);
    g = gifti(g);

    % convert vertices to original space
    vert_mat = double(([g.vertices ones(size(g.vertices,1),1)])');
    vert_mat = freeSurfer2T1*vert_mat;
    vert_mat(4,:) = [];
    vert_mat = vert_mat';
    g.vertices = vert_mat; clear vert_mat
       
    gifti_name = fullfile(dDir,subj,s_info.anat,['T1w_' fw_obj_hemi{kk} 'h_render.gii']);

    save(g,gifti_name,'Base64Binary')
end


%%
%% SVD on the odd responses:
%%

scan_nr = 3;

% load PPG responses
data_in = 'PPG';
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};

% nifti:
ni=niftiRead(fullfile(dDir,subj,scan,[scanName '.nii.gz']));

load(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponseT']),'t')

% load average of all odd heartbeats:
ppgTS=niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_odd.nii.gz']));

% load average of all odd heartbeats:
ppgTSeven=niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_even.nii.gz']));

% load coregistration matrix:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% scale the time-series matrix by the correlation

% Load the correlation with heartbeat (made with bbCorrelate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_cod' data_in '.nii.gz']);
ppgR = niftiRead(ppgRname); % COD between even and odd heartbeats

% Set maximum of ppgTS to 1 for each voxel
ppgTS.data = ppgTS.data ./ repmat(max(abs(ppgTS.data),[],4),[1,1,1,size(ppgTS.data,4)]);
ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
% Multiply by correlation size (absolute)
ppgTS.data = ppgTS.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTS.data,4)]));
ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));

% reshape to voxel X time
a = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);

% select times to include in SVD
if isequal(data_in,'PPG')
    a = a(:,t>=-.5 & t<=2);
    t_sel = t(t>=-.5 & t<=2);
elseif isequal(data_in,'RESP')
    a = a(:,t>=-.5 & t<=4);
    t_sel = t(t>=-.5 & t<=4);
end
meanTS = mean(a,2);
a = a-repmat(meanTS,1,size(a,2)); % subtract the mean
[u,s,v]=svd(a','econ');
s=diag(s);

%%%% test for the sign of the 2nd component (also check for 1st???)
% in the 2nd component, the first peak should be negative
[~,pm_i]=findpeaks(double(u(:,2)),'minpeakdistance',10);
[~,nm_i]=findpeaks(-double(u(:,2)),'minpeakdistance',10);
first_peak = min([pm_i; nm_i]);
if u(first_peak,2)<0
    % do nothing
elseif u(first_peak,2)>0
    % reverse sign pc and weights
    u(:,2)=-u(:,2); v(:,2) = -v(:,2);
end

% get to cumulative explained variance:
var_explained=cumsum(s.^2) / sum(s.^2)


%% Put SVD output in a structure
out = [];
for k=1:2
    out(k).weights = reshape(v(:,k),[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);
end

% %%%%% MODEL WITH 2 COMPONENTS:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';
svdResults.model = reshape(pred,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3) size(ppgTS.data,4)]);

%%%%% MODEL ERROR
% train and test-sets
train_set = reshape(ppgTS.data,[prod(ppgTS.dim(1:3)) ppgTS.dim(4)]);
test_set = reshape(ppgTSeven.data,[prod(ppgTSeven.dim(1:3)) ppgTSeven.dim(4)]);
% test-retest error
test_train_error = sqrt(sum((test_set - train_set).^2,2));
% model error
test_model_error = sqrt(sum((test_set - pred).^2,2));
% relative RMS error:
rel_rms_error = test_model_error./test_train_error;

svdResults.error = reshape(rel_rms_error,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);



%% Get functional voxels to plot with rendering

Rthreshold = .8;

% Make a mask for areas with some level of BOLD signal:
% brainTH = quantile(reshape(ni.data(:,:,:,1),1,[]),.75); 
% brainMask = ni.data(:,:,:,1)>brainTH;
brainTH = quantile(reshape(mean(ni.data(:,:,:,5:end),4),1,[]),.75); 
brainMask = mean(ni.data(:,:,:,5:end),4)>brainTH;

% Use mask:
[ii,jj,kk]=ind2sub(size(ppgR.data),find(ppgR.data>=Rthreshold & brainMask>0));
% No mask:
% [ii,jj,kk]=ind2sub(size(ppgR.data),find(ppgR.data>=Rthreshold));
ijk_func=[ii jj kk];
clear ii jj kk % housekeeping

xyz_anat = mrAnatXformCoords(acpcXform, ijk_func);

%%

figure
brainHandle=bbRenderGifti(g);
% brainHandle.FaceAlpha=.5; % make the brain transparent
plot3(xyz_anat(:,1),xyz_anat(:,2),xyz_anat(:,3),'r.','MarkerSize',10)
bbViewLight(270,0)





