        

%%
%% Calculate slope, peaktime and FWHM and save 
%%

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
        
    % Get base name
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

    % Load coregistration matrix (to write output in T1 space):
    load([save_name_base '_AcpcXform_new.mat']);
    acpcXform = acpcXform_new; clear acpcXform_new
    
    % load pc1
    pc1Weight = niftiRead([save_name_base '_space-T1w_canoPc1Weights.nii.gz']);
    % load pc2
    pc2Weight = niftiRead([save_name_base '_space-T1w_canoPc2Weights.nii.gz']);

    % load r weights
    r_weight = niftiRead([save_name_base '_space-T1w_canoPc12R.nii.gz']);
    
    % name for rel rms error weights
    rel_rms_error = niftiRead([save_name_base '_space-T1w_canoPc12RelRMSE.nii.gz']);
    
    % load canonical pcs
    load(fullfile(dDir,'derivatives','brainbeat','group',['canonicalPC_leavout' sub_labels{ss}]),'pc1','pc2','t_svd')

    % Get functional for physioGet
    ni = niftiRead(fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']));
    ppgT = load([save_name_base '_PPGtrigResponseT.mat'],'t');
    % Select times that were included in SVD: we want -0.5 to 1.5 heartbeat cycle
    physio      = physioCreate('nifti',ni);
    ppg_cycle   = 1./physioGet(physio,'PPGrate'); % length of 1 heartbeat cycle in secs
    t_sel       = ppgT.t(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle);

    % resample canonical pc1 to original timing:
    t_hr = linspace(min(t_sel),max(t_sel),128); % temp variable for resampling to t_sel
    y1 = interp1(t_hr,pc1,t_sel); 
    y2 = interp1(t_hr,pc2,t_sel); clear t_hr
    
    % line up beta weights
    beta_weights = [pc1Weight.data(:) pc2Weight.data(:)];

    % model for each voxel
    %%%% we can just use y1 and y2 to make these calculations in terms of time
    %%%% instead of percentage of heartbeat
    model_v = beta_weights(:,1:2)*[pc1';pc2']; % use 2 pcs
    excl_these = sum(model_v,2);

    % if pc1<0
    this_slope = diff(model_v,[],2);
    [max_slope,m1_ind] = max(this_slope,[],2);
    [min_slope,m2_ind] = min(this_slope,[],2);

    max_slope_larger = max_slope>abs(min_slope);

    % get all slopes lined up
    these_slopes = max_slope; 
    these_slopes(~max_slope_larger) = min_slope(~max_slope_larger);
    % check whether negative slope matches negative PC1 and reverse
    % length(find(these_slopes<0 & pc1Weight.data(:)>0))

    % get all times of steepest slope indices lined up, this is not the right
    % measure of latency though...
    these_times = m1_ind;
    these_times(~max_slope_larger) = m2_ind(~max_slope_larger);
    these_times = t_svd(these_times); % actually get seconds
    these_times(excl_these==0) = NaN; % zero slope non-interpretable time

    % get curve peak time 
    % get time of minimum if trough, maximum if peak - this means that slope
    % rules over PC1 weight, but note that there are only few voxels where this
    % is the case
    peak_t = NaN(size(these_slopes));
    for kk = 1:size(model_v,1)
        clear t_ind
        if these_slopes(kk)<0 % trough
            [~,t_ind] = min(model_v(kk,t_svd<.5)); % one cycle
            peak_t(kk) = t_svd(t_ind);
        elseif these_slopes(kk)>0 % peak
            [~,t_ind] = max(model_v(kk,t_svd<.5)); % one cycle
            peak_t(kk) = t_svd(t_ind);
        end
    end
    peak_t_perc = round(100*peak_t); %ranging -50:+50

    % calculate FWHM - error, negative for arteries
    peak_FWHM = NaN(size(these_slopes));
    for kk = 1:size(model_v,1)
        clear t_ind
        if these_slopes(kk)<0 % trough/min
            [this_min,t_ind] = min(model_v(kk,t_svd<.5)); % one cycle
            [this_max] = max(model_v(kk,t_svd<.5)); % one cycle
            this_th = .5*(this_min+this_max);
            steps_back = find(model_v(kk,t_ind:-1:1)>this_th,1)-1; % walk back from min
            steps_forw = find(model_v(kk,t_ind:1:end)>this_th,1)-1; % walk forward from min
            peak_FWHM(kk) = 100*(t_svd(steps_back+steps_forw)-t_svd(1));
        elseif these_slopes(kk)>0 % peak/max
            [this_max,t_ind] = max(model_v(kk,t_svd<.5)); % one cycle
            [this_min] = min(model_v(kk,t_svd<.5)); % one cycle
            this_th = .5*(this_min+this_max);
            steps_back = find(model_v(kk,t_ind:-1:1)<this_th,1)-1; % walk back from min
            steps_forw = find(model_v(kk,t_ind:1:end)<this_th,1)-1; % walk forward from min
            peak_FWHM(kk) = 100*(t_svd(steps_back+steps_forw)-t_svd(1));
        end
    end


    %%% write these niftis
    % add beta weights in niftis in functional space
    ni1 = ni; % slope pos/neg and strength
    ni1.data = ni1.data(:,:,:,1);
    ni1.data(:) = these_slopes;
    ni1.qto_xyz = acpcXform;
    ni1.qto_ijk = inv(acpcXform);
    ni1.sto_xyz = acpcXform;
    ni1.sto_ijk = inv(acpcXform);

    % name for positive or negative slope
    pc1_newName = [save_name_base '_space-T1w_modelSlope.nii.gz'];
    niftiWrite(ni1,pc1_newName)

    ni2 = ni; % onset
    ni2.data = ni2.data(:,:,:,1);
    ni2.data(:) = peak_t_perc;
    ni2.qto_xyz = acpcXform;
    ni2.qto_ijk = inv(acpcXform);
    ni2.sto_xyz = acpcXform;
    ni2.sto_ijk = inv(acpcXform);

    % name for time: ranging -50:+50
    pc2_newName = [save_name_base '_space-T1w_modelOnset.nii.gz'];
    niftiWrite(ni2,pc2_newName)

    ni3 = ni; % FWHM
    ni3.data = ni3.data(:,:,:,1);
    ni3.data(:) = peak_FWHM;
    ni3.qto_xyz = acpcXform;
    ni3.qto_ijk = inv(acpcXform);
    ni3.sto_xyz = acpcXform;
    ni3.sto_ijk = inv(acpcXform);

    % FWHM in percengate
    pc3_newName = [save_name_base '_space-T1w_model_FWHM.nii.gz']; 
    niftiWrite(ni3,pc3_newName)
end