
% 1: We take the principle components, resampled to the heartrate for all
% subjects.
% 2: Do a PCA to see how they replicate across subjects.
% 3: Run through the data to test how well each voxel is predicted by a
% combination of these two components.
%
% DH & BW 2018, vistalab

%% Save canonical heartbeat responses

% subjects indices
s_nr = [2 3 4];
% scan nrs, here take the 25 or 20 degree flip angle
% scan_nr  = [1 1 2];
scan_nr  = [3 3 3];

% load first two principle components for these subjects and scans
all_pcs = zeros(length(s_nr),2,128);
for kk = 1:length(s_nr)
    load(['./local/s-' int2str(s_nr(kk)) '_scan-' int2str(scan_nr(kk)) 'pc12'],'y1','y2','t_hr')
    all_pcs(kk,1,:) = y1;
    all_pcs(kk,2,:) = y2;
end
% reshape into matrix size subject/scan*2 X time
all_pcs = reshape(all_pcs,size(all_pcs,1)*size(all_pcs,2),size(all_pcs,3))';

% do the pca on the first two pca's
[u,s,v] = svd(all_pcs);
% s = diag(s);

t_svd = linspace(-.5,1.5,128);

% plot them:
figure
temp = u*s;
% plot(u(:,1:3))
plot(t_svd,temp(:,1:3))
xlabel('time (heartbeat cycles)')

pc1  = temp(:,1);
pc2  = temp(:,2);
pc3  = temp(:,3);

% save them:
% save(['./local/allsubs_pc12'],'pc1','pc2','pc3')


%% load canonical heartbeat responses and run through data
clear all
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% Load canonical heartbeat responses:
load(['./local/allsubs_pc12'],'pc1','pc2','pc3')

% Load PPG responses:
s_nr = 4;
scan_nr = 1;

s_info = bb_subs(s_nr);
subj=s_info.subj;
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};
data_in = 'PPG';

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

% Functionals:
ni=niftiRead(fullfile(dDir,subj,scan,[scanName '.nii.gz']));

% Load average of all odd heartbeats:
ppgTS=niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_odd.nii.gz']));

% Load average of all even heartbeats:
ppgTSeven=niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_even.nii.gz']));

%%%% Scale the time-series matrix by the reliability
% Get odd/even corr/corr (made with bbCod/Correlate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_cod' data_in '.nii.gz']);
% ppgRname = fullfile(dDir,subj,scan,[scanName '_corr' data_in '.nii.gz']);
ppgR = niftiRead(ppgRname); % COD between even and odd heartbeats
% Set maximum of ppgTS to 1 for each voxel
ppgTS.data = ppgTS.data ./ repmat(max(abs(ppgTS.data),[],4),[1,1,1,size(ppgTS.data,4)]);
ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
% Multiply by correlation size (absolute)
ppgTS.data = ppgTS.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTS.data,4)]));
ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));

% Load timing of heartbeat triggered responses:
load(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponseT']),'t')

% Load coregistration matrix:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% Times that were included in SVD:
physio      = physioCreate('nifti',ni);
ppg_cycle   = 1./physioGet(physio,'PPGrate');
t_sel = t(t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);

% Get even heartbeat responses that were included in SVD:
a = reshape(ppgTS.data,[numel(ppgTSeven.data(:,:,:,1)) length(t)]);
a = a(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);

% Resample canonical pc1 and pc2 to original timing:
t_hr = linspace(min(t_sel),max(t_sel),128);
y1 = interp1(t_hr,pc1,t_sel);
y2 = interp1(t_hr,pc2,t_sel);
% y3 = interp1(t_hr,pc3,t_sel);

disp('get beta weights')
beta_weights = zeros(size(a,1),2);
for kk = 1:size(a,1) % voxels
    if mod(kk,10000)==0, disp(['voxel ' int2str(kk) ' of ' int2str(size(a,1))]),end
    [B] = regress(a(kk,:)',[y1;y2]');
    beta_weights(kk,:) = B;
end
disp('done')

disp('get cod')
r_weights = zeros(size(a,1),1);
for kk = 1:size(a,1) % voxels
    if mod(kk,10000)==0, disp(['voxel ' int2str(kk) ' of ' int2str(size(a,1))]),end
    model_v = beta_weights(kk,:)*[y1;y2];
    r_weights(kk) = calccod(model_v',a(kk,:)',1,0,0)./100;
end
r_weights(isnan(r_weights)) = 0; % zero out NaN when model prediction is zeros
disp('done')

%% Plot R/Betas back on brain
 
% Quick brain mask:
brain_vect = ni.data(:,:,:,4);
brain_vect = brain_vect(:);

% Select voxels with decent reliability:
select_voxels = brain_vect>10 & r_weights>.7;
 
figure,hold on
subplot(1,2,1)
plot(beta_weights(select_voxels,1),beta_weights(select_voxels,2),'ko')
% subplot(1,2,2)
% plot3(beta_weights(select_voxels,1),beta_weights(select_voxels,2),beta_weights(select_voxels,3),'ko')

% set_1 = select_voxels & beta_weights(:,1)>=0; % CSF
% set_2 = select_voxels & beta_weights(:,1)<0 & beta_weights(:,2)>0; % veins
% set_3 = select_voxels & beta_weights(:,1)<0 & beta_weights(:,2)<0; % arteries
% plot(beta_weights(set_1,1),beta_weights(set_1,2),'go')
% plot(beta_weights(set_2,1),beta_weights(set_2,2),'ro')
% plot(beta_weights(set_3,1),beta_weights(set_3,2),'bo')
% axis equal

sliceThisDim = 2;
if s_nr == 1
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [0 26 17]; 
elseif s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-1 50 -21]; % for figures 2
%     curPos = [-10 -20 -21]; % for figures 1
%     curPos = [-11 34 -71]; % Carotid
%     curPos = [-2 26 -63]; % Basilar
%     curPos = [1 26 -21]; % SliceThisDim 1 Anterior Cerebral Artery, used in example
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
%     curPos = [0,4,38];
%     curPos = [0,4,38]; % for figure set
    curPos = [1 26 -63]; % x = 1 SliceThisDim 1 for Anterior Cerebral Artery
elseif s_nr == 4
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0 4 35];%[x x 38] % x = 0 nicely captures posterior arteries/veins, x = -10, anterior middle artery
end

ni_r = ni;
ni_r.data = ni_r.data(:,:,:,1);
ni_r.data(:) = r_weights;
ni_r.data(r_weights<.3) = 0;
ni_r.data(r_weights<0) = 0;
% bbOverlayDotsAnat(ni_r,niAnatomy,acpcXform,sliceThisDim,imDims,curPos);

ni_r = ppgR;
ni_r.data(ni_r.data<0.2) = 0;
% bbOverlayDotsAnat(ni_r,niAnatomy,acpcXform,sliceThisDim,imDims,curPos);

niColor = [1 0 0;0 0 1; 0 1 0];
R_th = 0.5;
niIntensity = ni;
niIntensity.data = zeros(size(ppgR.data));
niIntensity.data(ni_r.data<R_th) = 0;
niIntensity.data(beta_weights(:,1)>0 & r_weights>=R_th) = 3;
niIntensity.data(beta_weights(:,1)<0 & beta_weights(:,2)>0 & r_weights>=R_th) = 1;
niIntensity.data(beta_weights(:,1)<0 & beta_weights(:,2)<0 & r_weights>=R_th) = 2;

figure
sliceThisDim = 1;
curPos = [-1 17 28];%[-1/-9 18 28] % zero
bbOverlayDotsAnat_PickColor(niIntensity,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,niColor);

%%
ni1 = ni;
ni1.data = ni1.data(:,:,:,1);
ni1.data(:) = beta_weights(:,1);

ni2 = ni;
ni2.data = ni2.data(:,:,:,1);
ni2.data(:) = beta_weights(:,2);
bbOverlayDotsAnat_FancyColorCircle(ni1,ni2,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.7);


