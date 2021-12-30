
% 1: Select voxels with high COD and PC1<1 | PC1>1 
% 2: Get multi echo curves for these voxels
%
% DH 2020

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

%%
%% WIP
%%

%% load canonical heartbeat responses and run through data

% Load canonical heartbeat responses:
load(['./local/allsubs_pc12'],'pc1','pc2','pc3')

% Load PPG responses:
s_nr        = 5;
scan_nr     = 1;
s_info      = bb_subs(s_nr);
subj        = s_info.subj;
scan        = s_info.scan{scan_nr};
scanName    = s_info.scanName{scan_nr};
data_in     = 'PPG';

% Load canonical PCs from other subjects:
load(['./local/canonicalPC_leavout' int2str(s_nr)],'pc1','pc2','pc3','t_svd')

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

% Functionals:
ni = niftiRead(fullfile(dDir,subj,scan,[scanName '.nii.gz']));

% Load average of all odd heartbeats:
ppgTS = niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_odd.nii.gz']));

% Load average of all even heartbeats:
ppgTSeven = niftiRead(fullfile(dDir,subj,scan,[scanName '_' data_in 'trigResponse_even.nii.gz']));

%%%% Scale the time-series matrix by the reliability
% Get odd/even corr/corr (made with bbCod/Correlate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_cod' data_in '.nii.gz']);
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

% get even heartbeat responses that were included in SVD:
a = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);
a = a(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);
test_set = reshape(ppgTSeven.data,[prod(ppgTSeven.dim(1:3)) ppgTSeven.dim(4)]);
test_set = test_set(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);

% resample canonical pc1 and pc2 to original timing:
t_hr = linspace(min(t_sel),max(t_sel),128);
y1 = interp1(t_hr,pc1,t_sel);
y2 = interp1(t_hr,pc2,t_sel);
y3 = interp1(t_hr,pc3,t_sel);

% get beta values on PC1 and PC2 model for every voxel in even response using regression
disp('get beta weights')
beta_weights = zeros(size(a,1),2);
for kk = 1:size(a,1) % voxels
    if mod(kk,10000)==0, disp(['voxel ' int2str(kk) ' of ' int2str(size(a,1))]),end
%     [B] = regress(a(kk,:)',[y1;y2;y3]');
    [B] = regress(a(kk,:)',[y1;y2]');
    beta_weights(kk,:) = B;
end
disp('done'), clear B

% test whether model can predict odd responses better than even responses
disp('get cod')
r_weights = zeros(size(a,1),1);
relRMS_weights = zeros(size(a,1),1);
% error even versus odd responses
test_train_error = sqrt(sum((test_set - a).^2,2));
% model
model_v = beta_weights(:,1:2)*[y1;y2]; % use 2 pcs
% model_v = beta_weights*[y1;y2;y3];
% model error
test_model_error = sqrt(sum((test_set - model_v).^2,2));
% relative RMS error:
rel_rms_error = test_model_error./test_train_error;

% r_weights: model fits well
for kk = 1:size(a,1) % voxels
    if mod(kk,10000)==0, disp(['voxel ' int2str(kk) ' of ' int2str(size(a,1))]),end
    r_weights(kk) = calccod(model_v(kk,:)',a(kk,:)',1,0,0)./100;
end
r_weights(isnan(r_weights)) = 0; % zero out NaN when model prediction is zeros
disp('done')

% beta1 (pc1) and beta2 (pc2) 
ni1 = ni; % pc 1
ni1.data = ni1.data(:,:,:,1);
ni1.data(:) = beta_weights(:,1);
ni2 = ni; % pc 2
ni2.data = ni2.data(:,:,:,1);
ni2.data(:) = beta_weights(:,2);


%%
%% Get functional voxels with reliable R

Rthreshold = .7;% .3

% SPM segmentation
niSPM = niftiRead(fullfile(dDir,subj,scan,[scanName '_spmSeg.nii.gz']));
brainMask = niSPM.data>0;

% Use mask based on ppgR (even-odd)
select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);
% Use mask based on model fit
% select_voxels = find(r_weights>=Rthreshold & brainMask(:)>0);
% No mask:
% select_voxels = find(ppgR.data>=Rthreshold);

% Get indiced of selected voxels
[ii,jj,kk] = ind2sub(size(ppgR.data),select_voxels);
ijk_func = [ii jj kk];
clear ii jj kk % housekeeping

% Get xyz coordinates of voxels 
xyz_anat = mrAnatXformCoords(acpcXform, ijk_func);

% PC1/PC2 weights for voxels to plot:
pc12_rel = [ni1.data(select_voxels) ni2.data(select_voxels)];

%% now load ME data 

% gte multi echo data
if s_nr==4
    scan_nrs = {[4 5],[6 7]};
elseif s_nr==5
    scan_nrs = {[4 5],[6 7]};
elseif s_nr==6
    scan_nrs = {[4 5],[6 7],[9 10]};
end

scn_me = 1;%:length(scan_nrs)
% first echo
scan1 = s_info.scan{scan_nrs{scn_me}(1)};
scanName1 = s_info.scanName{scan_nrs{scn_me}(1)};

% load coregistration matrix:
me_acpcxform = load(fullfile(dDir,subj,scan1,[scanName1 'AcpcXform_new.mat']));

% load time T:
t_me = load(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponseT']),'t');

% load average of all heartbeats:
fname_s0 = [scanName1 '_PPGtrigResponse_S0.nii.gz'];
fname_t2s = [scanName1 '_PPGtrigResponse_T2s.nii.gz'];
ni_s0 = niftiRead(fullfile(dDir,subj,scan1,fname_s0));
ni_t2s = niftiRead(fullfile(dDir,subj,scan1,fname_t2s));

% Get indices of voxels in multi echo data
ijk_me = round(mrAnatXformCoords(inv(me_acpcxform.acpcXform_new), xyz_anat));

% Get timeseries (S0 T2*) of these voxels 
ts_s0 = NaN(size(ijk_me,1),size(ni_s0.data,4)); % voxels X t
ts_t2s = NaN(size(ijk_me,1),size(ni_s0.data,4)); % voxels X t
for kk = 1:size(ijk_me,1)
    ts_s0(kk,:) = ni_s0.data(ijk_me(kk,1),ijk_me(kk,2),ijk_me(kk,3),:);
    ts_t2s(kk,:) = ni_t2s.data(ijk_me(kk,1),ijk_me(kk,2),ijk_me(kk,3),:);
end


% outlier vector
outlier_sig = zeros(size(ts_s0,1),1);
outlier_sig(var(ts_t2s,[],2)>.1 | var(ts_s0,[],2)>.1) = 1;


%%

pc1_th = 0;

model_resp = 5*model_v(select_voxels,:);

figure('Position',[0 0 300 500])
subplot(2,1,1),hold on
% plot(t_sel,mean(model_resp(pc12_rel(:,1)<0 & outlier_sig==0,:)),':','Color',[.5 .5 .5])
y = 100*ts_t2s(pc12_rel(:,1)<0 & outlier_sig==0,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 1],'EdgeColor',[.5 .5 1],'FaceAlpha',.5)
plot(t_me.t,mean(y),'b','LineWidth',2)

y = 100*ts_s0(pc12_rel(:,1)<0 & pc12_rel(:,2)>=0 & outlier_sig==0,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 .5],'EdgeColor',[.5 .5 .5],'FaceAlpha',.5)
plot(t_me.t,mean(y),'k','LineWidth',2)

ylim([-3 3])
ylabel('% modulation')
xlabel('time (s)')
title('Negative PC1 (blood), black = S0, red = T2*')

subplot(2,1,2),hold on
% plot(t_sel,mean(model_resp(pc12_rel(:,1)>0  & outlier_sig==0,:)),':','Color',[.5 .5 .5])
y = 100*ts_t2s(pc12_rel(:,1)>0 & outlier_sig==0,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 1],'EdgeColor',[.5 .5 1],'FaceAlpha',.5)
plot(t_me.t,mean(y),'b','LineWidth',2)

y = 100*ts_s0(pc12_rel(:,1)>0 & outlier_sig==0,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 .5],'EdgeColor',[.5 .5 .5],'FaceAlpha',.5)
plot(t_me.t,mean(y),'k','LineWidth',2)

ylabel('% modulation')
xlabel('time (s)')
title('Positive PC1 (CSF), black = S0, blue = T2*')
ylim([-3 3])

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','me',['s' int2str(s_nr) '_scan' int2str(scan_nr) '_me' int2str(scan_nrs{scn}(1)) '_S0_T2s']))
% print('-painters','-r300','-depsc',fullfile(dDir,'figures','me',['s' int2str(s_nr) '_scan' int2str(scan_nr) '_me' int2str(scan_nrs{scn}(1)) '_S0_T2s']))
% 

%% Try to split slower and faster 
%%

figure('Position',[0 0 400 600])
subplot(4,1,1),hold on
plot([0 0],[-3 3],'k')
these_trials = pc12_rel(:,1)<0 & pc12_rel(:,2)>0 & outlier_sig==0;
plot(t_sel,mean(model_resp(these_trials,:)),':','Color',[.5 .5 .5])
y = 100*ts_s0(these_trials,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 .5],'EdgeColor',[.5 .5 .5])
plot(t_me.t,mean(y),'k')

y = 100*ts_t2s(these_trials,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 1],'EdgeColor',[.5 .5 1],'FaceAlpha',.5)
plot(t_me.t,mean(y),'b')
ylim([-3 3])
ylabel('% modulation')

subplot(4,1,2),hold on
plot([0 0],[-3 3],'k')
these_trials = pc12_rel(:,1)<0 & pc12_rel(:,2)<0 & outlier_sig==0;
y = 100*ts_s0(these_trials,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 .5],'EdgeColor',[.5 .5 .5])
plot(t_sel,mean(model_resp(these_trials,:)),':','Color',[.5 .5 .5])
plot(t_me.t,mean(y),'k')

y = 100*ts_t2s(pc12_rel(:,1)<0 & pc12_rel(:,2)<0 & outlier_sig==0,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 1],'EdgeColor',[.5 .5 1],'FaceAlpha',.5)
plot(t_me.t,mean(y),'b')
ylabel('% modulation')
ylim([-3 3])

subplot(4,1,3),hold on
plot([0 0],[-3 3],'k')
these_trials = pc12_rel(:,1)>0 & pc12_rel(:,2)<0 & outlier_sig==0;
plot(t_sel,mean(model_resp(these_trials,:)),':','Color',[.5 .5 .5])
y = 100*ts_s0(these_trials,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 .5],'EdgeColor',[.5 .5 .5])
plot(t_me.t,mean(y),'k')

y = 100*ts_t2s(these_trials,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 1],'EdgeColor',[.5 .5 1],'FaceAlpha',.5)
plot(t_me.t,mean(y),'b')
ylabel('% modulation')
ylim([-3 3])

subplot(4,1,4),hold on
plot([0 0],[-3 3],'k')
these_trials = pc12_rel(:,1)>0 & pc12_rel(:,2)>0 & outlier_sig==0;
plot(t_sel,mean(model_resp(these_trials,:)),':','Color',[.5 .5 .5])
y = 100*ts_s0(these_trials,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 .5],'EdgeColor',[.5 .5 .5])
plot(t_me.t,mean(y),'k')

y = 100*ts_t2s(these_trials,:);
y_low = mean(y)-2*std(y)/sqrt(size(y,1)); y_up = mean(y)+2*std(y)/sqrt(size(y,1));
fill([t_me.t t_me.t(end:-1:1)],[y_up y_low(end:-1:1)],[.5 .5 1],'EdgeColor',[.5 .5 1],'FaceAlpha',.5)
plot(t_me.t,mean(y),'b')
ylabel('% modulation')
ylim([-3 3])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'figures','me',['s' int2str(s_nr) '_scan' int2str(scan_nr) '_me' int2str(scan_nrs{scn_me}(1)) '_S0_T2s_earlylate']))
print('-painters','-r300','-depsc',fullfile(dDir,'figures','me',['s' int2str(s_nr) '_scan' int2str(scan_nr) '_me' int2str(scan_nrs{scn_me}(1)) '_S0_T2s_earlylate']))


%%
%% tests imagesc all S0 and T2* curves

figure
subplot(2,4,1)
imagesc(100*ts_s0(pc12_rel(:,1)<0 & pc12_rel(:,2)<0,:),[-3 3])
title('S0 pc1<0 pc2<0')
subplot(2,4,2)
imagesc(100*ts_s0(pc12_rel(:,1)<0 & pc12_rel(:,2)>0,:),[-3 3])
title('S0 pc1<0 pc2>0')

subplot(2,4,3)
imagesc(100*ts_t2s(pc12_rel(:,1)<0 & pc12_rel(:,2)<0,:),[-3 3])
subplot(2,4,4)
imagesc(100*ts_t2s(pc12_rel(:,1)<0 & pc12_rel(:,2)>0,:),[-3 3])

subplot(2,4,5)
imagesc(100*ts_s0(pc12_rel(:,1)>0 & pc12_rel(:,2)<0,:),[-3 3])
subplot(2,4,6)
imagesc(100*ts_s0(pc12_rel(:,1)>0 & pc12_rel(:,2)>0,:),[-3 3])

subplot(2,4,7)
imagesc(100*ts_t2s(pc12_rel(:,1)>0 & pc12_rel(:,2)<0,:),[-3 3])
subplot(2,4,8)
imagesc(100*ts_t2s(pc12_rel(:,1)>0 & pc12_rel(:,2)>0,:),[-3 3])

%%
figure
subplot(2,2,1)
imagesc(100*ts_s0(pc12_rel(:,1)<0,:),[-3 3])

subplot(2,2,2)
imagesc(100*ts_t2s(pc12_rel(:,1)<0,:),[-3 3])

subplot(2,2,3)
imagesc(100*ts_s0(pc12_rel(:,1)>0,:),[-3 3])

subplot(2,2,4)
imagesc(100*ts_t2s(pc12_rel(:,1)>0,:),[-3 3])
