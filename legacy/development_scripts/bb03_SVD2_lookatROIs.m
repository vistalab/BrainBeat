clear all
close all

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
dDir = '/biac4/wandell/data/BrainBeat/data';

% chdir(dDir)

%% The T2* data are here.  

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

s_nr = 2;
s_info = bb_subs(s_nr);
subj=s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

% % Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii.gz']));
% % load Veno rotation matrix:
% xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));


%%
%% now do an SVD on the odd responses:
%%

% load PPG responses
scan_nr = 2;
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};

% nifti:
ni=niftiRead(fullfile(dDir,subj,scan,[scanName]));

load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT']),'t')

% load average of all odd heartbeats:
ppgTS=niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_odd']));

% load average of all odd heartbeats:
ppgTSeven=niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_even']));

% load coregistration matrix:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))

% Load the correlation with heartbeat (made with bbCorrelate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_corrPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

% scale the time-series matrix by the correlation
% divide by max %%%% maybe think about z-scoring instead
ppgTS.data = ppgTS.data ./ repmat(max(abs(ppgTS.data),[],4),[1,1,1,size(ppgTS.data,4)]);
ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
% multiply by correlation size (absolute)
ppgTS.data = ppgTS.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTS.data,4)]));
ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));

% reshape for SVD
a = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);

%%%%% TEST LOW PASS FILTER START
% srate = 20;
% band = 2;
% Rp   = 1; Rs = 20; % third order Butterworth
% high_p =  band(1)*2/srate;
% delta = 0.001*2/srate;
% high_s = min(1-delta,high_p+0.1);
% [n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
% [bf_b,bf_a] = butter(n_band,wn_band,'low');
% for k=1:size(a,1)
%     a(k,:) = filtfilt(bf_b,bf_a,double(a(k,:)));
% end
%%%%% TEST END

a = a(:,t>=-.5 & t<=2);
t_sel = t(t>=-.5 & t<=2);
a = a-repmat(mean(a,2),1,size(a,2)); % subtract the mean
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

%%
%% put output in structures:
% put 2 components weights in a matrix
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

% save first component weight:
% ni_save = ni;
% ni_save.data = out(1).weights;
% ni_save.fname = fullfile(dDir,subj,scan,[scanName '_pc1_weights']);
% niftiWrite(ni_save,[ni_save.fname])

%% 
%% get both in same space
%%

if s_nr==2
    bb = [-90 -120 -120; 90 130 90]; 
elseif s_nr==3
    bb = [-90 -120 -100; 90 130 110]; 
end

mmPerVoxel = [1 1 1];
interpMethod = [0 0 0 0 0 0]; %nearest neighbour
% interpMethod = [1 1 1 0 0 0]; %trilinear

%%%% first funx to mm Anat space
img2reslice = ni.data(:,:,:,1);
xform = inv(acpcXform_new);
[newImg] = mrAnatResliceSpm(img2reslice, xform,   ...
                          bb,mmPerVoxel,interpMethod,...
                          false);

%%%% TS to mm Anat space
img2reslice = ppgTS.data;
xform = inv(acpcXform_new);
[newImgTS] = mrAnatResliceSpm(img2reslice, xform,   ...
                          bb,mmPerVoxel,interpMethod,...
                          false);

%%%% weigths to mm Anat space
interpMethod = [1 1 1 0 0 0]; %trilinear
img2reslice = out(1).weights;
xform = inv(acpcXform_new);
[newImgPC1] = mrAnatResliceSpm(img2reslice, xform,   ...
                          bb,mmPerVoxel,interpMethod,...
                          false);
img2reslice = out(2).weights;
xform = inv(acpcXform_new);
[newImgPC2] = mrAnatResliceSpm(img2reslice, xform,   ...
                          bb,mmPerVoxel,interpMethod,...
                          false);

interpMethod = [0 0 0 0 0 0]; %nearest neighbour
%%%% anat to mm Anat space
img2reslice = niAnatomy.data(:,:,:);
xform = inv(niAnatomy.qto_xyz);
[newImgAnat,xform] = mrAnatResliceSpm(img2reslice, xform,   ...
                          bb,mmPerVoxel,interpMethod,...
                          false);

%%%% ventricles to mm Anat space
allVentricles = niftiRead([s_info.freesurferDir 'nii/ventricles.nii.gz']);
img2reslice = allVentricles.data; 
xform = inv(allVentricles.qto_xyz);
[newImgVentricles] = mrAnatResliceSpm(img2reslice, xform,   ...
                          bb,mmPerVoxel,interpMethod,...
                          false);
    
%% get all rois PC and TS:

roisToSegmentNames = {...
    'L_lat_ventr'           %1
    'L_inferior_lat_ventr'  %2
    '3rd_ventr'             %3
    '4th_ventr'             %4
    'CSF'                   %5
    'L_choroid_plexus'      %6
    'R_lat_ventr'           %7
    'R_inferior_lat_ventr'  %8
    'R_choroid_plexus'      %9
    '5th_ventricle'};       %10

svdResults.roi = [];
tsVals = reshape(newImgTS,[numel(newImgTS(:,:,:,1)) length(t)]);

for k=1:length(roisToSegmentNames)
    roiName=[s_info.freesurferDir 'nii/' roisToSegmentNames{k} '.nii.gz'];
    currentVentricles = niftiRead(roiName);
    img2reslice = currentVentricles.data; 
    xform = inv(currentVentricles.qto_xyz);
    [newImgCVentricles] = mrAnatResliceSpm(img2reslice, xform,   ...
                              bb,mmPerVoxel,interpMethod,...
                              false);

    pc1Vals = newImgPC1(newImgCVentricles==1);
    pc2Vals = newImgPC2(newImgCVentricles==1);
    svdResults.roi(k).PCcurves = pc1Vals*u(:,1)' + pc2Vals*u(:,2)';
    
    svdResults.roi(k).TScurves = tsVals(newImgCVentricles==1,:);
    clear currentVentricles newImgCVentricles 
end
clear tsVals pc1Vals pc2Vals

%% plot curves all ROIS

roiColors = {[1 0 0],[0 1 1],[0 1 0],[0 0 1],[0 0 0],[1 0 1],[1 0 0],[0 1 1],[1 0 1],[0 1 .5]};
roiPlot = [1 7 3 4];

figure
for k=roiPlot
    subplot(2,1,1),hold on
    plot(t,median(svdResults.roi(k).PCcurves),'Color',roiColors{k})
end
ylabel('median model')
xlabel('time(s)')
legend(roisToSegmentNames(roiPlot))
for k=roiPlot
    subplot(2,1,2),hold on
    plot(t,median(svdResults.roi(k).TScurves),'Color',roiColors{k})
end
ylabel('median time series')
xlabel('time(s)')
legend(roisToSegmentNames(roiPlot))
% figure
% for k=[1 7 3 4 2 8 5]%1:length(roisToSegmentNames)
%     subplot(2,1,1),hold on
%     plot(t,mean(abs(svdResults.roi(k).PCcurves)),'Color',roiColors{k})
%     subplot(2,1,2),hold on
%     plot(t,mean(abs(svdResults.roi(k).TScurves)),'Color',roiColors{k})
% end

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',['./figures/ventricles/smooth_svd_' subj '_' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_pc12_curve'])

%% plot some ROIs over space

% region of interest to use
% roiName = ventricles; % all ventricles together
roisToSegmentNames = {...
    'L_lat_ventr'           %1
    'L_inferior_lat_ventr'  %2
    '3rd_ventr'             %3
    '4th_ventr'             %4
    'CSF'                   %5
    'L_choroid_plexus'      %6
    'R_lat_ventr'           %7
    'R_inferior_lat_ventr'  %8
    'R_choroid_plexus'      %9
    '5th_ventricle'};       %10
roisUse = [1 7 3 4];
newImgUse = zeros(size(newImgAnat));
for k = roisUse;
    roiName=[s_info.freesurferDir 'nii/' roisToSegmentNames{k} '.nii.gz']; 
    currentVentricles = niftiRead(roiName);
    img2reslice = currentVentricles.data; 
    xform = inv(currentVentricles.qto_xyz);
    [newImg] = mrAnatResliceSpm(img2reslice, xform,   ...
                              bb,mmPerVoxel,interpMethod,...
                              false);
    newImgUse(newImg>0) = 1;
    clear newImg
end

scale2max = .015;
figure,hold on
for x = 1:size(newImgAnat,1)
    for y = 1:size(newImgAnat,2)
        for z = 1:size(newImgAnat,3)
            if newImgUse(x,y,z)>0
            % now get the RGB values to plot from the weights
            val_plot = [newImgPC1(x,y,z) newImgPC2(x,y,z)];
            % scale to prespecified maximum:
            val_plot = val_plot./scale2max;
            % remove larger values to allow different max scaling
            val_plot(abs(val_plot)>1)=sign(val_plot(abs(val_plot)>1))*1;
            % get RGB
            data_colors_rgb = bbData2Colors(val_plot);
            plot3(x,y,z,'.','Color',data_colors_rgb)
            end
        end
    end
end
clear newImgUse
% view(80,40)
view(110,40)
axis equal
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/ventricles/smooth_svd_' subj '_' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_pc12_spatial'])

%%
%% plot clusters
%%

s_vect = sqrt(v(:,1).^2 + v(:,2).^2);

% quick brain mask:
brain_vect = ni.data(:,:,:,4);
brain_vect = brain_vect(:);

% correlation mask:
ppgR_vect = ppgR.data(:);

% select_voxels = s_vect>0.01 & brain_vect>10 & abs(ppgR_vect)>.3;
select_voxels = brain_vect>10 & ppgR_vect.^2>.7;

% get colors:
colors_plot = v(:,1:2);
colors_plot = colors_plot / max(colors_plot(:));
colors_plot = bbData2Colors(colors_plot);

figure
plot(v(select_voxels,1),v(select_voxels,2),'ko')
axis equal

v_plot = v(select_voxels,:);
colors_plot = colors_plot(select_voxels,:);
figure,hold on
for k=1:length(v_plot)
    plot(v_plot(k,1),v_plot(k,2),'o','Color',colors_plot(k,:))
end
axis equal


%%
%% reliability
%%

% model:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';

train_set = reshape(ppgTS.data,[prod(ppgTS.dim(1:3)) ppgTS.dim(4)]);
test_set = reshape(ppgTSeven.data,[prod(ppgTSeven.dim(1:3)) ppgTSeven.dim(4)]);

% error
test_train_error = sqrt(sum((test_set - train_set).^2,2));

% model error
test_model_error = sqrt(sum((test_set - pred).^2,2));

% relative RMS error:
rel_rms_error = test_model_error./test_train_error;

figure,hist(rel_rms_error,50)

%%
%% plot PC shapes and colors from data
%%

s_vect = sqrt(v(:,1).^2 + v(:,2).^2);

% model prediction for all voxels:
pred = [u(:,1:2)*diag(s(1:2))*v(:,1:2)']';

% quick brain mask:
brain_vect = ni.data(:,:,:,4);
brain_vect = brain_vect(:);

% correlation mask:
ppgR_vect = ppgR.data(:);

% select_voxels = s_vect>0.01 & brain_vect>10 & abs(ppgR_vect)>.3;
select_voxels = brain_vect>10 & ppgR_vect.^2>.5;

% get colors:
colors_plot = v(:,1:2);
colors_plot = colors_plot / max(abs(colors_plot(:)));
colors_plot = bbData2Colors(colors_plot);

figure
plot(v(select_voxels,1),v(select_voxels,2),'ko')
axis equal

v_plot = v(select_voxels,:);
voxel_nr = find(select_voxels>0);
colors_plot = colors_plot(select_voxels,:);

t_select = (t>-0.2 & t<1);

figure('Position',[0 0 700 700]),hold on
plot([-0.02 0.02],[0 0],'k')
plot([0 0],[-0.02 0.02],'k')
for k=1:length(v_plot) 
%     plot(v_plot(k,1),v_plot(k,2),'.','Color',colors_plot(k,:))

    plot(v_plot(k,1)+t(t_select)/2000,v_plot(k,2)+pred(voxel_nr(k),t_select)/2000,'Color',colors_plot(k,:))

end
axis equal

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/test_pca/svd_' subj '_pc12_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_plotcolorshape90'])

figure('Position',[0 0 300 300]),hold on
for k=1:length(v_plot) 
    plot(t(t_select),pred(voxel_nr(k),t_select),'Color',colors_plot(k,:))
end
axis tight
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',['./figures/test_pca/svd_' subj '_pc12_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_plotcolorshapeflat'])

%%
%% What do the colors mean?
%%

figure('Position',[0 0 700 700]),hold on
% example data and figure
data_in = [0 0;0 1;1 0;1 1;0 -1;-1 0;-1 -1;-1 1;1 -1];
data_in = [data_in;.2 * data_in; .4 * data_in; .6 * data_in; .8 * data_in];
data_colors_rgb = bbData2Colors(data_in);

for k=1:length(data_in)
%     plot(data_in(k,1),data_in(k,2),'k.','MarkerSize',40)
%     plot(data_in(k,1),data_in(k,2),'.','Color',data_colors_rgb(k,:),'MarkerSize',35)
    pred_thisval = data_in(k,1)*u(:,1) + data_in(k,2)*u(:,2);
    plot(data_in(k,1)+t(t_select)/3,data_in(k,2),'k','LineWidth',1)
    plot(data_in(k,1)+t(t_select)/3,data_in(k,2)+pred_thisval(t_select)/3,'Color',data_colors_rgb(k,:),...
        'LineWidth',2)
end

axis equal

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',['./figures/test_pca/svd_' subj '_pc12_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_exampleshapes'])

%%
figure('Position',[0 0 700 700]),hold on
% example data and figure
data_in = [0 0;0 1;1 0;1 1;0 -1;-1 0;-1 -1;-1 1;1 -1];
data_in = [data_in;.2 * data_in; .4 * data_in; .6 * data_in; .8 * data_in];
data_colors_rgb = bbData2Colors(data_in);

for k=1:length(data_in)
%     plot(data_in(k,1),data_in(k,2),'k.','MarkerSize',40)
%     plot(data_in(k,1),data_in(k,2),'.','Color',data_colors_rgb(k,:),'MarkerSize',35)
    pred_thisval = data_in(k,1)*u(:,1) + data_in(k,2)*u(:,2);
    plot(t(t_select),pred_thisval(t_select),'Color',data_colors_rgb(k,:),...
        'LineWidth',2)
end

axis equal

set(gcf,'PaperPositionMode','auto')
%     print('-painters','-r300','-dpng',fullfile(dDir,subj,scan,[scanName '_BBcurves_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))
print('-painters','-r300','-dpng',['./figures/test_pca/svd_' subj '_pc12_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_exampleshapes'])


