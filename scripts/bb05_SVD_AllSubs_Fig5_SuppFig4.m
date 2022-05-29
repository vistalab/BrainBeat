
% 
% This script is part of:
%%%%% Measuring brain beats: cardiac-aligned fast fMRI signals %%%%%
% Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell
%
%
% 0: We did the SVD in bb04_SVD.m for all individual subjects, here we
%   test across subjects.
% 1: We take the individual subject principle components, resampled to the
% heartrate for all subjects.
% 2: Do a PCA to see how they replicate across subjects.
% 3: Run through the data to test how well each voxel is predicted by a
% combination of these two components.
% 4: Save the beta weights for the canonical PCs in a nifti in T1 space
%
% DH & BW 2018, vistalab

clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%% Save canonical heartbeat responses (across all subjects)

sub_labels = {'1','2','3','4','5'}; 
ses_labels = {'1','1','1','1','1'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1]};

all_pcs = zeros(length(sub_labels),2,128);
all_varExplained = zeros(length(sub_labels),20);

for ss = 1:length(sub_labels) % subjects/ses/acq

    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);
    
    % Get PPG triggered curves
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);

    % load first two principle components for these subjects and scans
    load([save_name_base '_pc12'],'y1','y2','y3','t_hr','var_explained','all_pred_acc')
    % var explained is in training set, all_pred_acc is in testing set
    
    all_pcs(ss,1,:) = y1;
    all_pcs(ss,2,:) = y2;
    all_pcs(ss,3,:) = y3;

    all_varExplained(ss,:) = all_pred_acc(1:20); % just look at explained var from first 20 
end

t_svd = linspace(-.5,1.5,128);

% reshape into matrix size subject/scan*2 X time
all_pcs_temp = reshape(all_pcs,size(all_pcs,1)*size(all_pcs,2),size(all_pcs,3))';

% calculate the canonical PCs on the first two PCs
[u,s,v] = svd(all_pcs_temp);
temp = u*s;
pc1  = temp(:,1);
pc2  = temp(:,2);

% plot individual subject PCs:
figure('Position',[0 0 300 250])
subplot(2,2,1),hold on
plot(t_svd,squeeze(all_pcs(:,1,:)),'-','LineWidth',1)
plot(t_svd,squeeze(all_pcs(:,2,:)),'--','LineWidth',1)
title('PCs for 5 subjects')
xlim([t_svd(1) t_svd(end)])
set(gca,'XTick',[0 1],'YTick',[-0.2 0 0.2]),ylim([-.4 .5])

% plot individual subject PCs explained variance in test set
subplot(2,2,3),hold on
plot(all_varExplained','LineWidth',1)
xlim([0 20]),ylim([0 1]),set(gca,'XTick',[0:2:20],'YTick',[0 .25 .5 .75 1])

% plot canonical PCs
subplot(2,2,2),hold on
plot(t_svd,pc1,'k','LineWidth',1)
plot(t_svd,pc2,'k--','LineWidth',1)
xlim([t_svd(1) t_svd(end)])
xlabel('time (heartbeat cycles)')
title('canonical PCs')
set(gca,'XTick',[0 1],'YTick',[-0.5 0 0.5]),ylim([-.9 .9])

% fft of canonical PCs
L = length(temp(:,1)); % length
Fs = 1./mean(diff(t_svd)); % sampling frequency
f = Fs * (0:(L/2))/L; % frequency
p1 = abs(fft(temp(:,1:2))/L);
p1 = p1(1:floor(L/2+1),:); % first half
p1(2:end-1,:) = 2*p1(2:end-1,:); % 2 X first half except DC

subplot(2,2,4),hold on
plot(f,p1(:,1),'k','LineWidth',1)
plot(f,p1(:,2),'k--','LineWidth',1)
xlim([0 7])
ylabel('|P(f)|')
xlabel('frequency (Hz)')

set(gcf,'PaperPositionMode','auto')
print('-r300','-depsc2',fullfile(dDir,'derivatives','figures','Figure5AB_canonicalPC_S1-5'))
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures','Figure5AB_canonicalPC_S1-5'))

% save them:
save(fullfile(dDir,'derivatives','brainbeat','group','allsubs_pc12'),'pc1','pc2')

%% Figure 5 C and E and Supplemental Figure S4A Combinations of PC1 and PC2 
 
load(fullfile(dDir,'derivatives','brainbeat','group','allsubs_pc12'),'pc1','pc2')
pc1_plot = pc1(1:75);
pc2_plot = pc2(1:75);

figure('Position',[0 0 400 250])
subplot(2,2,1), plot(pc1_plot,'k','LineWidth',2), axis tight, axis off
title('canonical PCs')
subplot(2,2,3), plot(pc2_plot,'k','LineWidth',2), axis tight, axis off
[x,y] = meshgrid(-1:.4:1,-1:.4:1);
subplot(1,2,2),hold on
for kk = 1:size(x,1)
    for ll = 1:size(x,2)
        plot(x(kk,ll):.3/74:x(kk,ll)+.3, y(kk,ll)+.3*(x(kk,ll)*pc1_plot + y(kk,ll)*pc2_plot),'k','LineWidth',1)
    end
end
axis square, axis tight, title('predicted responses')
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures','Figure5CE_modelpc12'))
print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','figures','Figure5CE_modelpc12'))


figure('Position',[0 0 400 250])
subplot(2,2,1), plot(pc1_plot,'k','LineWidth',2), axis tight, axis off
title('canonical PCs')
subplot(2,2,3), plot(pc2_plot,'k','LineWidth',2), axis tight, axis off
[x,y] = meshgrid(-1:.4:1,-1:.4:1);
subplot(1,2,2),hold on
for kk = 1:size(x,1)
    for ll = 1:size(x,2)
        data_colors_rgb = bbData2Colors([x(kk,ll) y(kk,ll)]); % fancy 2D color map
        plot(x(kk,ll):.3/74:x(kk,ll)+.3, y(kk,ll)+.3*(x(kk,ll)*pc1_plot + y(kk,ll)*pc2_plot),...
            'Color',data_colors_rgb,'LineWidth',2)
    end
end
axis square, axis tight, title('predicted responses')
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures','SuppFigureS4A_modelpc12colors'))
print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','figures','SuppFigureS4A_modelpc12colors'))


%% Save canonical heartbeat responses (across N-1 subjects)

% run across the subjects and leave out subject kkk
for kk = 1:length(sub_labels)
    % leave our subj kk
    this_pcset = all_pcs(setdiff(1:length(sub_labels),kk),:,:);
    % reshape into matrix size subject/scan*2 X time
    Nmin1_pcs = reshape(this_pcset,size(this_pcset,1)*size(this_pcset,2),size(this_pcset,3))';

    % do the svd on pc1 and pc2
    [u,s,v] = svd(Nmin1_pcs);
    % s = diag(s);
    temp = u*s;
    pc1  = temp(:,1);
    pc2  = temp(:,2);
    
    % save them:
    save(fullfile(dDir,'derivatives','brainbeat','group',['canonicalPC_leavout' sub_labels{kk}]),'pc1','pc2','t_svd')

end

clear pc1 pc2 u s v temp Nmin1_pcs this_pcset

%% load canonical heartbeat responses and run through data
% now that we have canonical PCs, we can run these through the testing data
% and calculate relative RMSE.
% this cell:
% - saves beta weights on canonical PCs in a nifti in T1 space
% - saves COD (model-->data) and relative rms error (model VS test-retest) for canonical PC model

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
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']));
    % Load coregistration matrix (for the functionals)
    load([save_name_base '_AcpcXform_new.mat']);
    acpcXform = acpcXform_new; clear acpcXform_new

    % Get PPG triggered curves for training and testing set
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);
    ppgTSodd = niftiRead([save_name_base '_PPGtrigResponse_odd.nii.gz']); % ppg triggered time series
    ppgTSeven = niftiRead([save_name_base '_PPGtrigResponse_even.nii.gz']); % ppg triggered time series
    ppgT = load([save_name_base '_PPGtrigResponseT.mat'],'t');
    t = ppgT.t;

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
hist_Rrmse = zeros(length(sub_labels),2,length(0:.1:3));
out = [];

for ss = 1:length(sub_labels) % subjects/ses/acq

    rr = 1;% run_nr
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    run_nr = run_nrs{ss}(rr);
        
    % Get base name
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);

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
    [n,x] = hist(rel_rmse,0:.1:3);
    hist_Rrmse(ss,1,:) = n;
    hist_Rrmse(ss,2,:) = x;
    clear n x
    
    % fill outputs
    out(ss).pc1_th = pc1Weight.data(r_weight.data>0.7);
    out(ss).pc2_th = pc2Weight.data(r_weight.data>0.7);
    
    % get heartrate to have interpretable timing in seconds again
    ni = niftiRead(fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']));
    physio      = physioCreate('nifti',ni);
    out(ss).heartrate   = 60*physioGet(physio,'PPGrate'); %BpM
end


%% plot relative RMSE for all subjects
figure('Position',[0 0 100 100]),hold on
for ss = 1:5
    plot(squeeze(hist_Rrmse(ss,2,:)),squeeze(hist_Rrmse(ss,1,:)))
end
plot([1 1],[0 6000],'k')

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures','Figure5D_modelpc12_Rrmse'))
print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','figures','Figure5D_modelpc12_Rrmse'))


%% distribution of PC1/2 weights and plot with pc1/pc2 fancy latency color map
% circle size weighted by voxel density using hist3

% print heartrate
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
    disp(['sub ' int2str(ss) ' heartrate ' int2str(out(ss).heartrate) ' bpm'])
end

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures','SuppFigureS4B_modelpc12_weightsV2'))
print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','figures','SuppFigureS4B_modelpc12_weightsV2'))

