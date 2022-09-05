clear all
close all

%% Make a rendering 

% mris_convert lh.pial lh.pial.gii

%% Base data directory

% dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

s_nr = 6;
s_info = bb_subs(s_nr);
subj=s_info.subj;

hemi_load = 'r';

gifti_name = fullfile(dDir,subj,s_info.anat,['T1w_' hemi_load 'h_white_render.gii']);
% gifti_name = fullfile(dDir,subj,s_info.anat,['T1w_' hemi_load 'h_render.gii']);
if exist(gifti_name,'file')
    g = gifti(gifti_name);
else
    error('First make rendering for this subject')
end

%% Get the Flywheel obj file and make it a gifti:
if exist(gifti_name,'file')
    warning('Rendering already exists')
end
mri_orig = fullfile(dDir,subj,'freesurfer','mri','orig.mgz');
% fw_obj(1).name = fullfile(dDir,subj,'freesurfer','FlyWheelAnalyses','lh.pial.obj');
% fw_obj(2).name = fullfile(dDir,subj,'freesurfer','FlyWheelAnalyses','rh.pial.obj');
fw_obj(1).name = fullfile(dDir,subj,'freesurfer','FlyWheelAnalyses','lh.white.obj');
fw_obj(2).name = fullfile(dDir,subj,'freesurfer','FlyWheelAnalyses','rh.white.obj');
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
       
%     gifti_name = fullfile(dDir,subj,s_info.anat,['T1w_' fw_obj_hemi{kk} 'h_render.gii']);
    gifti_name = fullfile(dDir,subj,s_info.anat,['T1w_' fw_obj_hemi{kk} 'h_white_render.gii']);

    save(g,gifti_name,'Base64Binary')
end


%%
%% SVD on the odd responses:
%%

scan_nr = 1;

% load PPG responses
data_in = 'PPG';
scan = s_info.scan{scan_nr};
scanName = s_info.scanName{scan_nr};

% nifti:
ni = niftiRead(fullfile(dDir,subj,scan,[scanName '.nii.gz']));

load(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponseT']),'t')

% load average of all odd heartbeats:
ppgTS = niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_odd.nii.gz']));

% load average of all odd heartbeats:
ppgTSeven = niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_even.nii.gz']));

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
[u,s,v] = svd(a','econ');
s = diag(s);

%%%% test for the sign of the 2nd component (also check for 1st???)
% in the 2nd component, the first peak should be negative
[~,pm_i] = findpeaks(double(u(:,2)),'minpeakdistance',10);
[~,nm_i] = findpeaks(-double(u(:,2)),'minpeakdistance',10);
first_peak = min([pm_i; nm_i]);
if u(first_peak,2)<0
    % do nothing
elseif u(first_peak,2)>0
    % reverse sign pc and weights
    u(:,2) = -u(:,2); v(:,2) = -v(:,2);
end

% get to cumulative explained variance:
var_explained = cumsum(s.^2) / sum(s.^2)


%% Put SVD output in a structure
out = [];
for k = 1:2
    out(k).weights = reshape(v(:,k),[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);
end

%%%%% MODEL WITH 2 COMPONENTS:
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

Rthreshold = .6;

% Make a mask for areas with some level of BOLD signal in the average fMRI:
brainTH = quantile(reshape(mean(ni.data(:,:,:,5:end),4),1,[]),.75); 
brainMask = mean(ni.data(:,:,:,5:end),4)>brainTH;

% Use mask:
select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);
% No mask:
% select_voxels = find(ppgR.data>=Rthreshold);

% Get indiced of selected voxels
[ii,jj,kk] = ind2sub(size(ppgR.data),select_voxels);
ijk_func = [ii jj kk];
clear ii jj kk % housekeeping

% Get xyz coordinates of voxels 
xyz_anat = mrAnatXformCoords(acpcXform, ijk_func);

% PC1/PC2 weights for voxels to plot:
pc12_render = [out(1).weights(select_voxels) out(2).weights(select_voxels)];


% Select hemisphere
xyz_select = xyz_anat(:,1)>-10;
xx_plot = xyz_anat(xyz_select,1);
yy_plot = xyz_anat(xyz_select,2);
zz_plot = xyz_anat(xyz_select,3);

pc12_render_sel = pc12_render(xyz_select,:);

% Set maximum for dot colors:
maxPlot = .008;

% Make 2D colormap: one to vary color, the other varies intensity
cm = jet(250); cm = cm(26:225,:);
cm = cm(end:-1:1,:);
cm = cm+.4; cm(cm>1)=1;
gray_vect = .2*ones(200,3);
cm2D = zeros(100,size(cm,1),3);
for kk = 1:100
    cm2D(kk,:,:) = cm*kk/100 + gray_vect*(100-kk)/100;
end

% Get colors for selected voxels
intensity_plot = pc12_render_sel(:,1)./maxPlot;    
intensity_plot(intensity_plot>1) = 1;
intensity_plot(intensity_plot<-1) = -1;
color_plot = pc12_render_sel(:,2)./maxPlot;
color_plot(color_plot>1) = 1;
color_plot(color_plot<-1) = -1;

%% Plot negative PC1 (veins/arteries)
figure
brainHandle = bbRenderGifti(g); hold on
% brainHandle.FaceAlpha = .5; % Make the brain transparent

for kk = 1:size(intensity_plot,1)
    if intensity_plot(kk) < 0
        c_use = squeeze(cm2D(ceil(-intensity_plot(kk)*99+1),ceil(-color_plot(kk)*99.5)+100,:));
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',c_use)
    end
end

title(['R>' num2str(Rthreshold,3)])

% Just plot all points in red:
% plot3(xx_plot,yy_plot,zz_plot,'r.','MarkerSize',10)
bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/render/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1neg_rh_lat'])

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/render/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1neg_rh_med'])


% Plot positive PC1 (CSF)
figure
brainHandle=bbRenderGifti(g); hold on

for kk = 1:size(intensity_plot,1)
    if intensity_plot(kk) > 0
        c_use = squeeze(cm2D(ceil(intensity_plot(kk)*99+1),ceil(color_plot(kk)*99.5)+100,:));
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',c_use)
    end
end

title(['R>' num2str(Rthreshold,3)])

% Just plot all points in red:
% plot3(xx_plot,yy_plot,zz_plot,'r.','MarkerSize',10)
bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/render/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1pos_rh_lat'])

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/render/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1pos_rh_med'])

%% plain renders
figure
brainHandle=bbRenderGifti(g); hold on

bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/render/sub-' int2str(s_nr) '_rh_lat_pial'])

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/render/sub-' int2str(s_nr) '_rh_med_pial'])
