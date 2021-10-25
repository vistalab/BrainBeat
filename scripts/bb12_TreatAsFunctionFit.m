
% treating PC1 as a function does not quite work, steepness, width etc...
% better approach to describe with 2 pcs and then get accurate zero
% crossings, peak time and width from PCs

% resample canonical pc1 to original timing:
t_hr = linspace(min(t_sel),max(t_sel),128); % temp variable for resampling to t_sel
y1 = interp1(t_hr,pc1,t_sel); 
y2 = interp1(t_hr,pc2,t_sel); clear t_hr

%%% now we do our new thing where we vary Gain and Delta time
theseDelta = t_sel(t_sel>-0.25*ppg_cycle & t_sel<.25*ppg_cycle);
fitT = t_sel(t_sel>-0.2*ppg_cycle & t_sel<1.2*ppg_cycle);% length to fit
predY = zeros(length(fitT),length(theseDelta)); % time x deltas

for kk = 1:length(theseDelta)
    t_shift = t_sel+(theseDelta(kk));
    predY(:,kk) = y1(t_shift>-0.2*ppg_cycle & t_shift<1.2*ppg_cycle);
end
clear t_shift

fitMe = train_set(651,t_sel>=-0.2 & t_sel<=1.2*ppg_cycle);

figure
subplot(2,2,1),hold on
plot(fitT,fitMe,'r','LineWidth',2)
xlim([-0.5 2])

subplot(2,2,3),hold on
plot(fitT,predY)
hold on
plot(t_sel,y1,'k','LineWidth',2)
xlim([-0.5 2])

subplot(1,2,2),hold on

plot(theseDelta,predY'*fitMe')
xlabel('delta t')

% plot best fitting curve back on
subplot(2,2,1),hold on
[m_val,d_ind] = max(abs(predY'*fitMe'));

if m_val>0
    plot(fitT,predY(:,d_ind),'k:')
elseif m_val<0
    plot(fitT,-predY(:,d_ind),'k:')
end
        
    

%% Load canonical heartbeats across n-1 subjects and fit time and direction

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

for ss = 1%:length(sub_labels) % subjects/ses/acq

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

    % Load coregistration matrix (to save outputs in T1 space):
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
    ppg_cycle   = 1./physioGet(physio,'PPGrate'); % length of 1 heartbeat cycle in secs
    train_set   = train_set(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);
    test_set    = test_set(:,t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);
    t_sel       = t(t>=(0-(.5*ppg_cycle)) & t<=1.5*ppg_cycle);

    % resample canonical pc1 to original timing:
    t_hr = linspace(min(t_sel),max(t_sel),128); % temp variable for resampling to t_sel
    y1 = interp1(t_hr,pc1,t_sel); 
    y2 = interp1(t_hr,pc2,t_sel); clear t_hr

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

%     % name for pc1
%     pc1_newName = [save_name_base '_space-T1w_canoPc1Weights.nii.gz'];
%     niftiWrite(ni1,pc1_newName)
%     
%     % name for pc2
%     pc2_newName = [save_name_base '_space-T1w_canoPc2Weights.nii.gz'];
%     niftiWrite(ni2,pc2_newName)
% 
%     % name for r weights
%     r_newName = [save_name_base '_space-T1w_canoPc12R.nii.gz'];
%     niftiWrite(ni3,r_newName)
%     
%     % name for rel rms error weights
%     relRmse_newName = [save_name_base '_space-T1w_canoPc12RelRMSE.nii.gz'];
%     niftiWrite(ni4,relRmse_newName)
% 
%     % zero model
%     ni4.data(:) = zero_voxels;
%     zero_model = [save_name_base '_space-T1w_canoPC12ZeroModel.nii.gz'];
%     niftiWrite(ni4,zero_model)

    clear ni1 ni2 ni3 ni4 ni5
end
