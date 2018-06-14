clear all
close all

%% Make a rendering 

% mris_convert lh.pial lh.pial.gii

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

s_nr = 5;
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
% fw_obj(1).name = fullfile(dDir,subj,'freesurfer','FlywheelAnalyses','lh.pial.obj');
% fw_obj(2).name = fullfile(dDir,subj,'freesurfer','FlywheelAnalyses','rh.pial.obj');
fw_obj(1).name = fullfile(dDir,subj,'freesurfer','FlywheelAnalyses','lh.white.obj');
fw_obj(2).name = fullfile(dDir,subj,'freesurfer','FlywheelAnalyses','rh.white.obj');
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
svdResults = [];

for k=1:2
    out(k).weights = reshape(v(:,k),[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);
end

%%%%% Model prediction with 2 components:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';
svdResults.model = reshape(pred,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3) size(ppgTS.data,4)]);

%%%%% Model error with 2 components:
% train and test-sets
train_set = reshape(ppgTS.data,[prod(ppgTS.dim(1:3)) ppgTS.dim(4)]);
test_set = reshape(ppgTSeven.data,[prod(ppgTSeven.dim(1:3)) ppgTSeven.dim(4)]);
% test-retest error
test_train_error = sqrt(sum((test_set - train_set).^2,2));
% model error
test_model_error = sqrt(sum((test_set - pred).^2,2));
% relative RMS error:
rel_rms_error = test_model_error./test_train_error;
% put in the structure:
svdResults.error = reshape(rel_rms_error,[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);


%% Get functional voxels to plot with rendering

Rthreshold = .5;

% SPM segmentation
niSPM = niftiRead(fullfile(dDir,subj,scan,[scanName '_spmSeg.nii.gz']));
brainMask = niSPM.data>0;

% Use mask:
select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);
% No mask:
% select_voxels = find(ppgR.data>=Rthreshold);

% Get indiced of selected voxels
[ii,jj,kk]=ind2sub(size(ppgR.data),select_voxels);
ijk_func=[ii jj kk];
clear ii jj kk % housekeeping

% Get xyz coordinates of voxels 
xyz_anat = mrAnatXformCoords(acpcXform, ijk_func);

% PC1/PC2 weights for voxels to plot:
pc12_render = [out(1).weights(select_voxels) out(2).weights(select_voxels)];

% model for voxels to plot:
a = reshape(svdResults.model,[size(svdResults.model,1)* size(svdResults.model,2)* size(svdResults.model,3), size(svdResults.model,4)]);
model_render = a(select_voxels,:);
clear a

% Select hemisphere
if s_nr == 1 || s_nr == 2 || s_nr == 3 
    xyz_select = xyz_anat(:,1)>-10;
elseif s_nr == 4
    xyz_select = xyz_anat(:,1)>-15;
elseif s_nr == 5
    xyz_select = xyz_anat(:,1)>-4;
else
    disp('add subject number for x-hemisphere threshold')
end
xx_plot = xyz_anat(xyz_select,1);
yy_plot = xyz_anat(xyz_select,2);
zz_plot = xyz_anat(xyz_select,3);

pc12_render_sel = pc12_render(xyz_select,:);
model_render_sel = model_render(xyz_select,:); 

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
    if intensity_plot(kk)<0
        c_use = squeeze(cm2D(ceil(-intensity_plot(kk)*99+1),ceil(-color_plot(kk)*99.5)+100,:));
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',c_use)
    end
end

title(['R>' num2str(Rthreshold,3)])

% Just plot all points in red:
% plot3(xx_plot,yy_plot,zz_plot,'r.','MarkerSize',10)
% bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/render/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1neg_rh_lat'])

bbViewLight(270,0)
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/render/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1neg_rh_med'])

% Plot positive PC1 (CSF)
figure
brainHandle = bbRenderGifti(g); hold on

for kk = 1:size(intensity_plot,1)
    if intensity_plot(kk)>0
        c_use = squeeze(cm2D(ceil(intensity_plot(kk)*99+1),ceil(color_plot(kk)*99.5)+100,:));
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',c_use)
    end
end

title(['R>' num2str(Rthreshold,3)])

% Just plot all points in red:
% plot3(xx_plot,yy_plot,zz_plot,'r.','MarkerSize',10)
% bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/render/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1pos_rh_lat'])

bbViewLight(270,0)
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/render/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1pos_rh_med'])

%% render 1 timepoint

figure
brainHandle = bbRenderGifti(g); hold on
% brainHandle.FaceAlpha = .5; % Make the brain transparent

cm = jet(100);

% move to indices 1:100 for colormap:
model_render_plot = model_render_sel - min(model_render_sel(:)); % set minimum to 0
model_render_plot = round(99*(model_render_plot./max(model_render_plot(:))))+1;

for tt = 1%:length(t)
    for kk = 1:size(model_render_sel,1)
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',cm(model_render_plot(kk,tt),:))
    end
end
title(['R>' num2str(Rthreshold,3)])

% Just plot all points in red:
% plot3(xx_plot,yy_plot,zz_plot,'r.','MarkerSize',10)
bbViewLight(90,0)

%% Upsample t and timeseries

% move to indices 1:100 for colormap:
model_render_plot = model_render_sel - min(model_render_sel(:)); % set minimum to 0
model_render_plot = round(99*(model_render_plot./max(model_render_plot(:))))+1;

upFactor = 8;
tUp = interp(t,upFactor);
modelUp = zeros(size(model_render_plot,1),length(tUp));

for kk = 1:size(model_render_plot,1)
    modelUp(kk,:) = interp(double(model_render_plot(kk,:)),upFactor);
        %imgSliceUp(kx,ky,:) = smooth(imgSliceUp(kx,ky,:),20); %maybe smooth
end
% correct for outliers after interpolation
modelUp(modelUp<1) = 1;
modelUp(modelUp>100) = 100;
    
% Get 1 heartbeat cycle:
physio      = physioCreate('nifti',ni);
ppg_cycle   = 1./physioGet(physio,'PPGrate');
t_start = find(tUp>=-ppg_cycle/2,1);
t_end = find(tUp>ppg_cycle/2,1);

%% movie

videoName =  [dDir 'movies/Render/model_sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_viewRH_2_med'];

% video 1:
% cm = hot(100);
% video 2:
cm = jet(100);

fid = figure('Position',[0 0 500 500]);

vidObj = VideoWriter(videoName,'MPEG-4'); %

open(vidObj); 

for tt = t_start:t_end % loop over time
    
    brainHandle = bbRenderGifti(g); hold on
    
    for kk = 1:size(model_render_sel,1)
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',cm(round(modelUp(kk,tt)),:))
    end
    
%     bbViewLight(90,0)
    bbViewLight(270,0)
    
    
    title(['t = ' num2str(tUp(tt),3)])
    
    % Let the number of frames to write depend on the timing, such that
    % peaks can be emphasized
    nr_frames = 1;
%     if t(kk)>-0.3 && t(kk)<0.3 
%         nr_frames = 20;
%     else
%         nr_frames = 5;
%     end
    
    % Write each frame to the file
    for m = 1:nr_frames % write X frames: decides speed
        writeVideo(vidObj,getframe(fid));
    end
    
    clf
    
end

close(vidObj);

%% negative pulses
videoName =  [dDir 'movies/Render/model_sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_viewRH_med_Neg_3'];

% video 1:
% cm = hot(100);
% video 2:
% cm = jet(100);
% video 3:
cm = hot(42);
cm = [cm(1,:)+zeros(8,3); cm; cm(end,:)+zeros(50,3)];

figure('Position',[0 0 500 500]);
brainHandle1 = bbRenderGifti(g); hold on
bbViewLight(270,0)
for kk = 1:size(model_render_sel,1)
    plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'k.','MarkerSize',30)
end
ax = gca;

fid = figure('Position',[0 0 800 500]);
vidObj = VideoWriter(videoName,'MPEG-4'); %
open(vidObj); 
for tt = t_start:t_end % loop over time
    
    brainHandle1 = bbRenderGifti(g); 
    bbViewLight(270,0)
    set(gca,'Xlim',ax.XLim,'Ylim',ax.YLim,'Zlim',ax.ZLim), hold on
    title(['t = ' num2str(tUp(tt),3)])

    for kk = 1:size(model_render_sel,1)
        if intensity_plot(kk) < 0 && round(modelUp(kk,tt))<35
            plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),...
                '.','MarkerSize',30,'Color',cm(round(modelUp(kk,tt)),:))
        end            
    end

    
    % Let the number of frames to write depend on the timing, such that
    % peaks can be emphasized
    nr_frames = 5;
%     if t(kk)>-0.3 && t(kk)<0.3 
%         nr_frames = 20;
%     else
%         nr_frames = 5;
%     end
    
    % Write each frame to the file
    for m = 1:nr_frames % write X frames: decides speed
        writeVideo(vidObj,getframe(fid));
    end
    
    cla
    
end

close(vidObj);

%% positive pulses
videoName =  [dDir 'movies/Render/model_sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_viewRH_med_Pos_3'];

% video 1:
% cm = hot(100);
% video 2:
% cm = jet(100);
% video 3:
cm = jet(100);
cm = cm(end:-1:1,:); 

figure('Position',[0 0 500 500]);
brainHandle1 = bbRenderGifti(g); hold on
bbViewLight(270,0)
for kk = 1:size(model_render_sel,1)
    plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'k.','MarkerSize',30)
end
ax = gca;

fid = figure('Position',[0 0 800 500]);
vidObj = VideoWriter(videoName,'MPEG-4'); %
open(vidObj); 
for tt = t_start:t_end % loop over time
    
    brainHandle1 = bbRenderGifti(g); 
    bbViewLight(270,0)
    set(gca,'Xlim',ax.XLim,'Ylim',ax.YLim,'Zlim',ax.ZLim), hold on
    title(['t = ' num2str(tUp(tt),3)])

    for kk = 1:size(model_render_sel,1)
        if intensity_plot(kk)>0 && round(modelUp(kk,tt))>65
            plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),...
                '.','MarkerSize',30,'Color',cm(floor(modelUp(kk,tt)),:))
        end            
    end

    
    % Let the number of frames to write depend on the timing, such that
    % peaks can be emphasized
    nr_frames = 5;
%     if t(kk)>-0.3 && t(kk)<0.3 
%         nr_frames = 20;
%     else
%         nr_frames = 5;
%     end
    
    % Write each frame to the file
    for m = 1:nr_frames % write X frames: decides speed
        writeVideo(vidObj,getframe(fid));
    end
    
    cla
    
end

close(vidObj);


%% negative pulses Lateral
videoName =  [dDir 'movies/Render/model_sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_viewRH_lat_Neg_3'];

% video 1:
% cm = hot(100);
% video 2:
% cm = jet(100);
% video 3:
cm = hot(42);
cm = [cm(1,:)+zeros(8,3); cm; cm(end,:)+zeros(50,3)];

figure('Position',[0 0 500 500]);
brainHandle1 = bbRenderGifti(g); hold on
bbViewLight(90,0)
for kk = 1:size(model_render_sel,1)
    plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'k.','MarkerSize',30)
end
ax = gca;

fid = figure('Position',[0 0 800 500]);
vidObj = VideoWriter(videoName,'MPEG-4'); %
open(vidObj); 
for tt = t_start:t_end % loop over time
    
    brainHandle1 = bbRenderGifti(g); 
    bbViewLight(90,0)
    set(gca,'Xlim',ax.XLim,'Ylim',ax.YLim,'Zlim',ax.ZLim), hold on
    title(['t = ' num2str(tUp(tt),3)])

    for kk = 1:size(model_render_sel,1)
        if intensity_plot(kk) < 0 && round(modelUp(kk,tt))<35
            plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),...
                '.','MarkerSize',30,'Color',cm(round(modelUp(kk,tt)),:))
        end            
    end

    
    % Let the number of frames to write depend on the timing, such that
    % peaks can be emphasized
    nr_frames = 5;
%     if t(kk)>-0.3 && t(kk)<0.3 
%         nr_frames = 20;
%     else
%         nr_frames = 5;
%     end
    
    % Write each frame to the file
    for m = 1:nr_frames % write X frames: decides speed
        writeVideo(vidObj,getframe(fid));
    end
    
    cla
    
end

close(vidObj);


%% positive pulses lateral
videoName =  [dDir 'movies/Render/model_sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_viewRH_lat_Pos_3'];

% video 1:
% cm = hot(100);
% video 2:
% cm = jet(100);
% video 3:
cm = jet(100);
cm = cm(end:-1:1,:); 

figure('Position',[0 0 500 500]);
brainHandle1 = bbRenderGifti(g); hold on
bbViewLight(90,0)
for kk = 1:size(model_render_sel,1)
    plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'k.','MarkerSize',30)
end
ax = gca;

fid = figure('Position',[0 0 800 500]);
vidObj = VideoWriter(videoName,'MPEG-4'); %
open(vidObj); 
for tt = t_start:t_end % loop over time
    
    brainHandle1 = bbRenderGifti(g); 
    bbViewLight(90,0)
    set(gca,'Xlim',ax.XLim,'Ylim',ax.YLim,'Zlim',ax.ZLim), hold on
    title(['t = ' num2str(tUp(tt),3)])

    for kk = 1:size(model_render_sel,1)
        if intensity_plot(kk)>0 && round(modelUp(kk,tt))>65
            plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),...
                '.','MarkerSize',30,'Color',cm(floor(modelUp(kk,tt)),:))
        end            
    end

    
    % Let the number of frames to write depend on the timing, such that
    % peaks can be emphasized
    nr_frames = 5;
%     if t(kk)>-0.3 && t(kk)<0.3 
%         nr_frames = 20;
%     else
%         nr_frames = 5;
%     end
    
    % Write each frame to the file
    for m = 1:nr_frames % write X frames: decides speed
        writeVideo(vidObj,getframe(fid));
    end
    
    cla
    
end

close(vidObj);

%% All pulses
videoName =  [dDir 'movies/Render/model_sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_viewRH_med_All_3B'];

% video 1:
% cm = hot(100);
% video 2:
cm = jet(100);

figure('Position',[0 0 500 500]);
brainHandle1 = bbRenderGifti(g); hold on
bbViewLight(270,0)
for kk = 1:size(model_render_sel,1)
    plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'k.','MarkerSize',20)
end
ax = gca;

fid = figure('Position',[0 0 500 500]);
vidObj = VideoWriter(videoName,'MPEG-4'); %
open(vidObj); 

for tt = t_start:t_end % loop over time
    
    brainHandle1 = bbRenderGifti(g); 
    bbViewLight(270,0)
    set(gca,'Xlim',ax.XLim,'Ylim',ax.YLim,'Zlim',ax.ZLim), hold on
    title(['t = ' num2str(tUp(tt),3)])

    for kk = 1:size(model_render_sel,1)
%         if intensity_plot(kk)<0 && round(modelUp(kk,tt))<35 
%             plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),...
%                 '.','MarkerSize',20,'Color',cm(round(modelUp(kk,tt)),:))
%         elseif intensity_plot(kk)>0 && round(modelUp(kk,tt))>65
%             plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),...
%                 '.','MarkerSize',20,'Color',cm(round(modelUp(kk,tt)),:))
%         end  
        if round(modelUp(kk,tt))<25 || round(modelUp(kk,tt))>75
            plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),...
                '.','MarkerSize',20,'Color',cm(round(modelUp(kk,tt)),:))
        end  
    end
    
    % Let the number of frames to write depend on the timing, such that
    % peaks can be emphasized
    nr_frames = 5;
%     if t(kk)>-0.3 && t(kk)<0.3 
%         nr_frames = 20;
%     else
%         nr_frames = 5;
%     end
    
    % Write each frame to the file
    for m = 1:nr_frames % write X frames: decides speed
        writeVideo(vidObj,getframe(fid));
    end
    
    cla
    
end

close(vidObj);

%% Plot some shapes and colorbars

figure('Position',[0 0 100 200])
% negative responses
negResp = modelUp(intensity_plot<0,:);
cm = hot(42);
cm = [cm(1,:)+zeros(8,3); cm; cm(end,:)+zeros(50,3)];
% plot(tUp,negResp(1,:))

subplot(2,1,1),hold on
nn = 277; % response number to plot
plot(tUp,negResp(nn,:),'Color',[.5 .5 .5],'LineWidth',2)
for kk = 1:length(tUp)
    if negResp(nn,kk)<35
        plot(tUp(kk),negResp(nn,kk),'.','Color',cm(round(negResp(nn,kk)),:),'MarkerSize',10)
    end
end
% nn = 853; % response number to plot
% plot(tUp,negResp(nn,:),'Color',[.5 .5 .5],'LineWidth',2)
% for kk = 1:length(tUp)
%     if negResp(nn,kk)<35
%         plot(tUp(kk),negResp(nn,kk),'.','Color',cm(round(negResp(nn,kk)),:),'MarkerSize',10)
%     end
% end
xlim([tUp(t_start) tUp(t_end)])

% positive responses
posResp = modelUp(intensity_plot>0,:);
cm = jet(100);
cm = cm(end:-1:1,:); 
axis off

subplot(2,1,2),hold on
nn = 594; % response number to plot
plot(tUp,posResp(nn,:),'Color',[.5 .5 .5],'LineWidth',2)
for kk = 1:length(tUp)
    if posResp(nn,kk)>65
        plot(tUp(kk),posResp(nn,kk),'.','Color',cm(round(posResp(nn,kk)),:),'MarkerSize',10)
    end
end
% nn = 340; % response number to plot
% plot(tUp,posResp(nn,:),'Color',[.5 .5 .5],'LineWidth',2)
% for kk = 1:length(tUp)
%     if posResp(nn,kk)>65
%         plot(tUp(kk),posResp(nn,kk),'.','Color',cm(round(posResp(nn,kk)),:),'MarkerSize',10)
%     end
% end
xlim([tUp(t_start) tUp(t_end)])
axis off
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-depsc',[dDir './movies/render/colormap'])
print('-painters','-r300','-dpng',[dDir './movies/render/colormap'])

