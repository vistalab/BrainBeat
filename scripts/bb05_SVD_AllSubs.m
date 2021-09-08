
% 1: We take the principle components, resampled to the heartrate for all
% subjects.
% 2: Do a PCA to see how they replicate across subjects.
% 3: Run through the data to test how well each voxel is predicted by a
% combination of these two components.
% 4: Save the beta weights for the canonical PCs in a nifti in T1 space
%
% DH & BW 2018, vistalab

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

%% Save canonical heartbeat responses (across all subjects)

sub_labels = {'1','2','3','4','5'}; 
ses_labels = {'1','1','1','1','1'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1]};

all_pcs = zeros(length(sub_labels),2,128);

for ss = 1:length(sub_labels) % subjects/ses/acq

    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);
    
    % Get PPG triggered curves
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

    % load first two principle components for these subjects and scans
    load([save_name_base '_pc12'],'y1','y2','y3','t_hr','var_explained','all_pred_acc')

    all_pcs(ss,1,:) = y1;
    all_pcs(ss,2,:) = y2;
end

t_svd = linspace(-.5,1.5,128);

% plot them:
figure('Position',[0 0 300 250])
subplot(2,1,1),hold on
plot(t_svd,squeeze(all_pcs(:,1,:)),'Color',[0 0 .5],'LineWidth',1)
plot(t_svd,squeeze(all_pcs(:,2,:)),'Color',[1 .5 0],'LineWidth',1)
% plot(t_svd,squeeze(all_pcs(:,3,:)),'Color',[1 .5 1],'LineWidth',2)
title('heartbeat components for 6 subjects')
xlim([t_svd(1) t_svd(end)])
set(gca,'XTick',[0 1],'YTick',[-0.2 0 0.2])

% reshape into matrix size subject/scan*2 X time
all_pcs_temp = reshape(all_pcs,size(all_pcs,1)*size(all_pcs,2),size(all_pcs,3))';

% do the pca on the first two pca's
[u,s,v] = svd(all_pcs_temp);
% s = diag(s);
temp = u*s;
pc1  = temp(:,1);
pc2  = temp(:,2);

subplot(2,2,3),hold on
plot(t_svd,pc1,'Color',[0 0 .5],'LineWidth',2)
plot(t_svd,pc2,'Color',[1 .5 0],'LineWidth',2)
xlim([t_svd(1) t_svd(end)])
xlabel('time (heartbeat cycles)')
title('canonical heartbeat components')
set(gca,'XTick',[0 1],'YTick',[-0.5 0 0.5])

% fft of pc
L = length(temp(:,1)); % length
Fs = 1./mean(diff(t_svd)); % sampling frequency
f = Fs * (0:(L/2))/L; % frequency
p1 = abs(fft(temp(:,1:2))/L);
p1 = p1(1:floor(L/2+1),:); % first half
p1(2:end-1,:) = 2*p1(2:end-1,:); % 2 X first half except DC

subplot(2,2,4),hold on
plot(f,p1(:,1),'Color',[0 0 .5],'LineWidth',2)
plot(f,p1(:,2),'Color',[1 .5 0],'LineWidth',2)
xlim([0 7])
ylabel('|P(f)|')
xlabel('frequency (Hz)')

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','canonicalPC_S1-5'))
print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group','canonicalPC_S1-5'))

% save them:
save(fullfile(dDir,'derivatives','brainbeat','group','allsubs_pc12'),'pc1','pc2')

%% Combinations of PC1 and PC2 with color scale:

load(fullfile(dDir,'derivatives','brainbeat','group','allsubs_pc12'),'pc1','pc2')

pc1_plot = pc1(1:75);
pc2_plot = pc2(1:75);

figure('Position',[0 0 400 250])
subplot(2,2,1)
plot(pc1_plot,'k','LineWidth',2)
axis tight
axis off
subplot(2,2,3)
plot(pc2_plot,'k','LineWidth',2)
axis tight
axis off

% Make 2D fancy colormap
[x,y] = meshgrid(-1:.4:1,-1:.4:1);
data_in = [x(:) y(:)];
data_colors_rgb = bbData2Colors(data_in);

subplot(1,2,2),hold on
for kk = 1:size(data_in,1)
    x = data_in(kk,1);
    y = data_in(kk,2);
    plot([x:.3/74:x+.3],y+.3*(x*pc1_plot + y*pc2_plot),...
        'Color',data_colors_rgb(kk,:),...
        'LineWidth',2)
end
axis square
axis tight
% axis off
title('canonical PCs across 5 subjects')

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','modelpc12'))
print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group','modelpc12'))


%% Save canonical heartbeat responses (across N-1 subjects)

% run across the subjects and leave out subject kkk
for kk = 1:length(sub_labels)
    % leave our subj kk
    this_pcset = all_pcs(setdiff(1:length(sub_labels),kk),:,:);
    % reshape into matrix size subject/scan*2 X time
    Nmin1_pcs = reshape(this_pcset,size(this_pcset,1)*size(this_pcset,2),size(this_pcset,3))';

    % do the pca on the first two pca's
    [u,s,v] = svd(Nmin1_pcs);
    % s = diag(s);
    temp = u*s;
    pc1  = temp(:,1);
    pc2  = temp(:,2);
%     pc3  = temp(:,3);
    
    % save them:
    save(fullfile(dDir,'derivatives','brainbeat','group',['canonicalPC_leavout' sub_labels{kk}]),'pc1','pc2','t_svd')

end

clear pc1 pc2 u s v temp Nmin1_pcs this_pcset

%% load canonical heartbeat responses and run through data
% saves beta weights on canonical PCs in a nifti in T1 space
% saves COD (model-->data) and rms error (model VS test-retest) for canonical PC model

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

for ss = 1:length(sub_labels) % subjects/ses/acq

    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);

    % Load canonical PCs from other subjects - these are calculated on odd
    % heartbeats
    load(fullfile(dDir,'derivatives','brainbeat','group',['canonicalPC_leavout' sub_labels{ss}]),'pc1','pc2','t_svd')

    % Get functional for physioGet
    ni = niftiRead(fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']));

    % Get PPG triggered curves
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

    ppgTSodd = niftiRead([save_name_base '_PPGtrigResponse_odd.nii.gz']); % ppg triggered time series
    ppgTSeven = niftiRead([save_name_base '_PPGtrigResponse_even.nii.gz']); % ppg triggered time series
    ppgT = load([save_name_base '_PPGtrigResponseT.mat'],'t');
    t = ppgT.t;

    % Load coregistration matrix (for the functionals):
    load([save_name_base '_AcpcXform_new.mat']);
    acpcXform = acpcXform_new; clear acpcXform_new

    %%%% COD 
    ppgRname = [save_name_base '_codPPG.nii.gz'];
    ppgR = niftiRead(ppgRname); % correlation with PPG

    % Set maximum of ppgTS to 1 for each voxel
    ppgTSodd.data = ppgTSodd.data ./ repmat(max(abs(ppgTSodd.data),[],4),[1,1,1,size(ppgTSodd.data,4)]);
    ppgTSeven.data = ppgTSeven.data ./ repmat(max(abs(ppgTSeven.data),[],4),[1,1,1,size(ppgTSeven.data,4)]);
    % Multiply by COD size (absolute)
    ppgTSodd.data = ppgTSodd.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSodd.data,4)]));
    ppgTSeven.data = ppgTSeven.data .* abs(repmat(ppgR.data,[1,1,1,size(ppgTSeven.data,4)]));

    % Reshape to voxel X time:
    train_set = reshape(ppgTSodd.data,[numel(ppgTSodd.data(:,:,:,1)) length(t)]);
    % test data
    test_set = reshape(ppgTSeven.data,[numel(ppgTSeven.data(:,:,:,1)) length(t)]);

    % Select times that were included in SVD: we want -0.5 to 1.5 heartbeat cycle
    physio      = physioCreate('nifti',ni);
    ppg_cycle   = 1./physioGet(physio,'PPGrate');
    train_set   = train_set(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);
    test_set    = test_set(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);
    t_sel       = t(t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);

    % resample canonical pc1 and pc2 to original timing:
    t_hr = linspace(min(t_sel),max(t_sel),128);
    y1 = interp1(t_hr,pc1,t_sel);
    y2 = interp1(t_hr,pc2,t_sel);

    % get beta values on PC1 and PC2 model for every voxel in even response using regression
    disp('get beta weights')
    beta_weights = zeros(size(train_set,1),2);
    for kk = 1:size(train_set,1) % voxels
        if mod(kk,10000)==0, disp(['voxel ' int2str(kk) ' of ' int2str(size(train_set,1))]),end
        [B] = regress(train_set(kk,:)',[y1;y2]');
        beta_weights(kk,:) = B;
    end
    disp('done')

    % test whether model can predict odd responses better than even responses
    disp('get cod')
    r_weights = zeros(size(train_set,1),1);
    % error even versus odd responses
    test_train_error = sqrt(sum((test_set - train_set).^2,2));
    % model
    model_v = beta_weights(:,1:2)*[y1;y2]; % use 2 pcs
    % model error
    test_model_error = sqrt(sum((test_set - model_v).^2,2));
    % relative RMS error:
    rel_rms_error = test_model_error./test_train_error;

    for kk = 1:size(train_set,1) % voxels
        if mod(kk,10000)==0, disp(['voxel ' int2str(kk) ' of ' int2str(size(train_set,1))]),end
        r_weights(kk) = calccod(model_v(kk,:)',train_set(kk,:)',1,0,0)./100;
    end
    r_weights(isnan(r_weights)) = 0; % zero out NaN when model prediction is zeros
    disp('done')

    % get rid of zero voxels - zero model;
    zero_voxels = sum(test_set,2)==0;
    
    % save beta weights in a nifti in functional space 
    % bb08_MNI_PCA writes it to T1 and then to MNI space
    % this should be corrected, because we now save here in T1 space

    % add beta weights in niftis in functional space
    ni1 = ni; % PC1
    ni1.data = ni1.data(:,:,:,1);
    ni1.data(:) = beta_weights(:,1);
    ni1.qto_xyz = acpcXform;
    ni1.qto_ijk = inv(acpcXform);
    ni1.sto_xyz = acpcXform;
    ni1.sto_ijk = inv(acpcXform);

    ni2 = ni; % PC2
    ni2.data = ni2.data(:,:,:,1);
    ni2.data(:) = beta_weights(:,2);
    ni2.qto_xyz = acpcXform;
    ni2.qto_ijk = inv(acpcXform);
    ni2.sto_xyz = acpcXform;
    ni2.sto_ijk = inv(acpcXform);

    ni3 = ni; % r_weights
    ni3.data = ni3.data(:,:,:,1);
    ni3.data(:) = r_weights;
    ni3.qto_xyz = acpcXform;
    ni3.qto_ijk = inv(acpcXform);
    ni3.sto_xyz = acpcXform;
    ni3.sto_ijk = inv(acpcXform);

    ni4 = ni; % rel rms error
    ni4.data = ni4.data(:,:,:,1);
    ni4.data(:) = rel_rms_error;
    ni4.qto_xyz = acpcXform;
    ni4.qto_ijk = inv(acpcXform);
    ni4.sto_xyz = acpcXform;
    ni4.sto_ijk = inv(acpcXform);

    % name for pc1
    pc1_newName = [save_name_base '_space-T1w_canoPc1Weights.nii.gz'];
    niftiWrite(ni1,pc1_newName)
    
    % name for pc2
    pc2_newName = [save_name_base '_space-T1w_canoPc2Weights.nii.gz'];
    niftiWrite(ni2,pc2_newName)

    % name for r weights
    r_newName = [save_name_base '_space-T1w_canoPc12R.nii.gz'];
    niftiWrite(ni3,r_newName)
    
    % name for rel rms error weights
    relRmse_newName = [save_name_base '_space-T1w_canoPc12RelRMSE.nii.gz'];
    niftiWrite(ni4,relRmse_newName)

    % zero model
    ni4.data(:) = zero_voxels;
    zero_model = [save_name_base '_space-T1w_canoPC12ZeroModel.nii.gz'];
    niftiWrite(ni4,zero_model)

    clear ni1 ni2 ni3 ni4 ni5
end

%% Summary across subjects

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

modelBetterThanData = zeros(length(sub_labels),1);
out = [];

for ss = 1:length(sub_labels) % subjects/ses/acq

    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);

    % Get base name
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

    % load pc1
    pc1Weight = niftiRead([save_name_base '_space-T1w_canoPc1Weights.nii.gz']);
    % load pc2
    pc2Weight = niftiRead([save_name_base '_space-T1w_canoPc2Weights.nii.gz']);

    % load r weights
    r_weight = niftiRead([save_name_base '_space-T1w_canoPc12R.nii.gz']);
    
    % name for rel rms error weights
    rel_rms_error = niftiRead([save_name_base '_space-T1w_canoPc12RelRMSE.nii.gz']);
    
    % zero data and model
    zero_model = niftiRead([save_name_base '_space-T1w_canoPC12ZeroModel.nii.gz']);
    
    % voxels with model->data better than data->data
    rel_rmse = rel_rms_error.data(:);
    rel_rmse(zero_model.data==1) = [];
    modelBetterThanData(ss) = 100*length(find(rel_rmse<1))./length(rel_rmse);

    % fill outputs
    out(ss).pc1_th = pc1Weight.data(r_weight.data>0.7);
    out(ss).pc2_th = pc2Weight.data(r_weight.data>0.7);
end

%% plot for every subject in 1 color
cm = lines(length(sub_labels));
figure,hold on
for ss = 1:length(sub_labels) % subjects/ses/acq
    subplot(1,6,ss)
    p1 = scatter(out(ss).pc1_th,out(ss).pc2_th,5,cm(ss,:));
    set(p1,'MarkerFaceAlpha',.5)
    set(p1,'MarkerEdgeAlpha',.5)
end

%%
%% plot with pc1/pc2 fancy latency color map
% circle size weighted by voxel density using hist3

figure('Position',[0 0 800 200])
for ss = 1:length(sub_labels) % subjects/ses/acq
    X = [out(ss).pc1_th,out(ss).pc2_th];

    [n,c] = hist3(X,'Ctrs',{-1.5:0.2:1.5 -1.5:0.2:1.5});
    subplot(1,6,ss),hold on
    for kk = 1:size(n,1)
        for ll = 1:size(n,2)
            
            data_in = [c{1}(kk) c{2}(ll)];
            data_in(data_in>1) = 1;
            data_in(data_in<-1) = -1;
            data_colors_rgb = bbData2Colors(data_in);
            
            plot(c{1}(kk),c{2}(ll),'.','MarkerSize',1+40*(n(kk,ll)/max(n(:))),...
                'Color',data_colors_rgb)
            xlim([-1.6 1.6]),ylim([-1.6 1.6])
            axis square
            
        end
    end
end
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','modelpc12_weightsV2'))
% print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group','modelpc12_weightsV2'))


%% plot 1 voxel to check model versus data
voxel_nr = find(ppgR.data(:)>.8,1);
figure('Position',[0 0 200 150]),hold on
plot(t_sel,model_v(voxel_nr,:),'k','LineWidth',2)
plot(t_sel,a(voxel_nr,:),'b:','LineWidth',2)
plot(t_sel,test_set(voxel_nr,:),'r:','LineWidth',2)

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/renderCanonSvd/subj' int2str(s_nr) '_scan' int2str(scan_nr) '_pred_train_test'])
% print('-painters','-r300','-depsc',[dDir './figures/renderCanonSvd/subj' int2str(s_nr) '_scan' int2str(scan_nr) '_pred_train_test'])

%%
%% FIGURES & RENDERINGS in individual subjects
%%
%% Plot R/Betas back on brain
%%
%%

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

ss = 6;
rr = 1;% run_nr
sub_label = sub_labels{ss};
ses_label = ses_labels{ss};
acq_label = acq_labels{ss};
run_nr = run_nrs{ss}(rr);

% Get base name of saved data
save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

% load pc1
pc1Weight = niftiRead([save_name_base '_space-T1w_canoPc1Weights.nii.gz']);
% load pc2
pc2Weight = niftiRead([save_name_base '_space-T1w_canoPc2Weights.nii.gz']);
% load r weights (how much explained by canonical PC1-2)
r_weight = niftiRead([save_name_base '_space-T1w_canoPc12R.nii.gz']);

% Get anatomy
t1w_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'anat',...
            ['sub-' sub_label '_ses-' ses_label '_T1w.nii.gz']);
niAnatomy = niftiRead(fullfile(dDir,t1w_BIDSname));

sliceThisDim = 1;
if ss == 1
    imDims = [-90 -120 -120; 90 130 90];
    curPos = [-4 26 17]; 
elseif ss == 2
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [-1 50 -21]; 
elseif ss == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [-2 26 -63]; 
elseif ss == 4
    imDims = [-90 -120 -50; 90 130 120];
    curPos = [0 4 35];
elseif ss == 5
    imDims = [-90 -120 -100; 90 130 120];
    curPos = [-2,18,38];
elseif ss == 6
    imDims = [-90 -120 -100; 90 130 120];
    curPos = [-4,18,38];
end

% plot beta1 (pc1) and beta2 (pc2) using fancy color circle
acpcXform = pc1Weight.qto_xyz;

% plot entire circle
bbOverlayDotsAnat_FancyColorCircle(pc1Weight,pc2Weight,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.7);
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_exampleSag']))
% print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_exampleSag']))
% 
% % now only select 7.30-1.30 on a clock
% % do this to test:
% % x = [1:-0.1:0 0:-0.1:-1 -1:0.1:0 0:0.1:1];
% % y = [0:0.1:1 1:-0.1:0 0:-0.1:-1 -1:0.1:0];
% % b = complex(x,y)
% % figure,plot(angle(b*(1-i)))
% 
% pc_complex = complex(pc1Weight.data,pc2Weight.data);
% pc_angle = angle(pc_complex*(1-i)); % multiply by (1-i) to rotatio 45 deg
% 
% % plot veins and arteries
% pc1_plot = pc1Weight;
% pc2_plot = pc2Weight;
% pc1_plot.data(pc_angle<0) = 0;
% pc2_plot.data(pc_angle<0) = 0;
% bbOverlayDotsAnat_FancyColorCircle(pc1_plot,pc2_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.7);
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_example_pc12veins']))
% print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_example_pc12veins']))
% 
% % plot csf
% pc1_plot = pc1Weight;
% pc2_plot = pc2Weight;
% pc1_plot.data(pc_angle>0) = 0;
% pc2_plot.data(pc_angle>0) = 0;
% bbOverlayDotsAnat_FancyColorCircle(pc1_plot,pc2_plot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.7);
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_example_pc12csf']))
% print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_example_pc12csf']))

%%
%% Get functional voxels to plot with rendering

Rthreshold = .5;

% get a segmentation to create a brain, CSF and vessels mask:
% niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
niSegm2 = niftiRead([save_name_base '_combineSegm.nii.gz']);
brainMask = niSegm2.data>0;

% load PPG-RCOD
ppgR = niftiRead([save_name_base '_codPPG.nii.gz']);

% Use mask:
select_voxels = find(ppgR.data>=Rthreshold & brainMask>0);

% Get indiced of selected voxels
[ii,jj,kk] = ind2sub(size(ppgR.data),select_voxels);
ijk_func = [ii jj kk];
clear ii jj kk % housekeeping

% Get xyz coordinates of voxels in single subject space
xyz_anat = mrAnatXformCoords(acpcXform, ijk_func);

% PC1/PC2 weights for voxels to plot:
pc12_render = [pc1Weight.data(select_voxels) pc2Weight.data(select_voxels)];

% Render right hemisphere
hemi_load = 'r';
save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

gifti_name = fullfile(dDir,'derivatives','surfaces',['sub-' sub_label],['ses-' ses_label],['T1w_' hemi_load 'h_white_render.gii']);
g = gifti(gifti_name);

% get coordinates within only this hemipshere
% xyz_select = xyz_anat(:,1)<10;
if ss==2 || ss==3
    x_plot = -15; % midline off in this subject, should have acpc oriented...
else
    x_plot = -10;
end
xyz_select = xyz_anat(:,1)>x_plot;
xx_plot = xyz_anat(xyz_select,1);
yy_plot = xyz_anat(xyz_select,2);
zz_plot = xyz_anat(xyz_select,3);

pc12_render_sel = pc12_render(xyz_select,:);

maxPlot = 0.5;
intensity_plot = pc12_render_sel./maxPlot;    
intensity_plot(intensity_plot>1) = 1;
intensity_plot(intensity_plot<-1) = -1;
data_colors_rgb = bbData2Colors([intensity_plot(:,1) intensity_plot(:,2)]);

figure
brainHandle = bbRenderGifti(g); hold on
for kk = 1:size(intensity_plot,1)
    plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',data_colors_rgb(kk,:))
end
% brainHandle.FaceAlpha = .5; % Make the brain transparent
title(['R>' num2str(Rthreshold,3)])

bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewLat']))

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewMed']))

% Render left hemisphere
hemi_load = 'l';
save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

gifti_name = fullfile(dDir,'derivatives','surfaces',['sub-' sub_label],['ses-' ses_label],['T1w_' hemi_load 'h_white_render.gii']);
g = gifti(gifti_name);

% get coordinates within only this hemipshere
if ss==2 || ss==3
    x_plot = 5; % midline off in this subject, should have acpc oriented...
else
    x_plot = 10;
end
xyz_select = xyz_anat(:,1)<x_plot;
xx_plot = xyz_anat(xyz_select,1);
yy_plot = xyz_anat(xyz_select,2);
zz_plot = xyz_anat(xyz_select,3);

pc12_render_sel = pc12_render(xyz_select,:);

maxPlot = 0.5;
intensity_plot = pc12_render_sel./maxPlot;    
intensity_plot(intensity_plot>1) = 1;
intensity_plot(intensity_plot<-1) = -1;
data_colors_rgb = bbData2Colors([intensity_plot(:,1) intensity_plot(:,2)]);

figure
brainHandle = bbRenderGifti(g); hold on
for kk = 1:size(intensity_plot,1)
    plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',data_colors_rgb(kk,:))
end
% brainHandle.FaceAlpha = .5; % Make the brain transparent
title(['R>' num2str(Rthreshold,3)])

bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewMed']))

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['subj' int2str(ss) '_run' int2str(rr) '_render' upper(hemi_load) '_viewLat']))



%%
%%
%%
%%
%% only plot upper left/lower righ quadrant
%%
pc_complex = complex(intensity_plot(:,1),intensity_plot(:,2));
pc_angle = angle(pc_complex*(1-1i)); % multiply by (1-i) to rotatio 45 deg
figure
brainHandle = bbRenderGifti(g); hold on
for kk = 1:size(intensity_plot,1)
    if pc_angle(kk)<0
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',data_colors_rgb(kk,:))
    end
end
title(['R>' num2str(Rthreshold,3)])


%%
%% plot with 2D red/blue colormap
%%

% Set maximum for dot colors:
maxPlot = .5;

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

bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','renderCanonSvd',['sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1neg_rh_lat']))

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','renderCanonSvd',['sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1neg_rh_med']))

%% Plot positive PC1 (CSF)
figure
brainHandle = bbRenderGifti(g); hold on

for kk = 1:size(intensity_plot,1)
    if intensity_plot(kk)>0
        c_use = squeeze(cm2D(ceil(intensity_plot(kk)*99+1),ceil(color_plot(kk)*99.5)+100,:));
        plot3(xx_plot(kk),yy_plot(kk),zz_plot(kk),'.','MarkerSize',20,'Color',c_use)
    end
end

title(['R>' num2str(Rthreshold,3)])

bbViewLight(90,0)
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','renderCanonSvd',['sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1pos_rh_lat']))

bbViewLight(270,0)
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','renderCanonSvd',['sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_PC1pos_rh_med']))

