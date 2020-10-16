function [response_matrix1,t,response_matrix_odd1,response_matrix_even1,response_matrix_std1,...
    response_matrix2,response_matrix_odd2,response_matrix_even2,response_matrix_std2] = ...
        bbResponse2physioME(ni1,ni2,slices,varargin)

% function to get brain response after the peak of a heartbeat
% meant to deal with multi echo (ME) data
%
% [response_matrix,t,response_matrix_odd,response_matrix_even] = bbResponse2physioME(ni,slices)
% 
% ni1: a nifti structure loaded by niftiRead of the first echo
% ni2: a nifti structure loaded by niftiRead of the second echo
%
% optional inputs:
% slices: slices for which to calculate response function, if no input, do
% all slices
%
% [response_matrix1,t,response_matrix_odd1,response_matrix_even1,response_matrix_std1,...
%     response_matrix2,response_matrix_odd2,response_matrix_even2,response_matrix_std2] = ...
%     bbResponse2physioME(ni,slices)
%
% DHermes, Copyright Vistasoft Team 2014

if ~exist('slices','var') % do whole brain
    slices = [1:size(ni1.data,3)];
elseif exist('slices','var') && isempty(slices)% do whole brain
    slices = [1:size(ni1.data,3)];
end

% get the nifti stuff we need from the first echo, same for echo 2:
timing = bbGet(ni1,'timing');
mux_f = bbGet(ni1,'super slices');
srate = 1/bbGet(ni1,'tr');

% get physio stuff we need:
physio     = physioCreate('nifti',ni1);
if isempty(varargin) % do PPG
    ppg_onsets = physioGet(physio,'ppg peaks');
elseif ~isempty(varargin) && isequal(varargin{2},'ppg') % do PPG
    ppg_onsets = physioGet(physio,'ppg peaks');
elseif ~isempty(varargin) && isequal(varargin{2},'resp') % do RESP
    ppg_onsets = physioGet(physio,'resp peaks');
end
    
% set epoch times
if isempty(varargin)
    epoch_pre = .5;%sec pre-onset
    epoch_post = 2;%sec post-onset
else
    epoch_pre = varargin{1}(1);%sec pre-onset
    epoch_post = varargin{1}(2);%sec post-onset
end
step_size = 1/srate/mux_f;% in s
srate_epochs = 1/step_size;
t = [-epoch_pre:step_size:epoch_post];%s timing for 1 epoch

% adding 1 time point in case there 
t_vox = min(timing(:)):step_size:(max(timing(:)+step_size)); % maximal timing accuracy for this scan session

% only use those ppg onsets that have an entire trial before and after
ppg_onsets = ppg_onsets((ppg_onsets-epoch_pre)>1); % get rid of early ones, and the first scans
ppg_onsets = ppg_onsets((ppg_onsets+epoch_post)<max(timing(:))); % get rid of late ones

% initiate output matrices to fill for s0 and T2*
response_matrix1 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_std1 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_odd1 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_even1 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));

response_matrix2 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_std2 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_odd2 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_even2 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));


for s = 1:length(slices)
    disp(['slice ' int2str(s) ' of ' int2str(length(slices))])
    sli = slices(s);
    % get the log for this slice, so we can solve the linear form
    d1 = log(squeeze(ni1.data(:,:,sli,1:end)));
    d2 = log(squeeze(ni2.data(:,:,sli,1:end)));

    % get echo times
    te1 = bbGet(ni1,'te');
    te2 = bbGet(ni2,'te');
    % calculate T2*
    rr = (d2 - d1)/(te2 - te1);
    t2s = -1./rr;
    
    % calculate S0
    lns0 = d1 - (d2 - d1)*(te1/(te2 - te1));
    s0 = exp(lns0);

    % plot to check for subject 3, scan 4-5
%     figure
%     subplot(2,1,1), hold on
%     plot(squeeze(d1(27,50,:)),'b')
%     plot(squeeze(d2(27,50,:)),'g')
%     plot(ppg_onsets*srate,8.5,'r.')    
%     subplot(2,1,2), hold on, ylabel('zscore')
%     plot(zscore(squeeze(s0(27,50,:))),'k')
%     plot(zscore(squeeze(t2s(27,50,:))),'m')
%     plot(ppg_onsets*srate,0,'r.')

    
    % demean and express in percent modulation 
    d1_norm = reshape(s0,[size(s0,1) * size(s0,2), size(s0,3)]);
    d2_norm = reshape(t2s,[size(t2s,1) * size(t2s,2), size(t2s,3)]);
    points_use = 4:size(d1_norm,2); % do not use the first couple of scans
    for kk = 1:size(d1_norm,1)
        % do not use first scans:
        x = points_use;
        
        % S0
        y = d1_norm(kk,points_use);       
        % detrend
        p = polyfit(x,y,1);    
        % percent modulation
        mean_factor = mean(d1_norm(kk,points_use)); % mean
        d1_norm(kk,:) = (d1_norm(kk,:) - (p(1)*[1:size(d1_norm,2)] + p(2)))./mean_factor;

        % T2*
        y = d2_norm(kk,points_use);       
        % detrend
        p = polyfit(x,y,1);    
        % percent modulation
        mean_factor = mean(d2_norm(kk,points_use)); % mean
        d2_norm(kk,:) = (d2_norm(kk,:) - (p(1)*[1:size(d2_norm,2)] + p(2)))./mean_factor;

    end
    d1 = reshape(d1_norm,[size(d1,1), size(d1,2), size(d1,3)]);
    d2 = reshape(d2_norm,[size(d2,1), size(d2,2), size(d2,3)]);
    clear d1_norm d2_norm

    % create a set of NaNs at the upsampled time rate 
    d1_up = single(NaN(size(d1,1),size(d1,2),length(t_vox))); % initiate with NaNs
    d2_up = single(NaN(size(d2,1),size(d2,2),length(t_vox))); % initiate with NaNs

    % get the timing for this slice (times where we have data)
    t_sli = timing(sli,:);

    % fill the data into these times:
    % get the indices in upsampled time at which data were collected
    d1_up(:,:,ismember(round(t_vox,3),round(t_sli,3))) = d1;
    d2_up(:,:,ismember(round(t_vox,3),round(t_sli,3))) = d2;

    % temporary response matrix for 1 slice: voxel X voxel X time X PPGonset
    temp_response_matrix1 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(t),length(ppg_onsets)));
    temp_response_matrix2 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(t),length(ppg_onsets)));
    % run through all ppg onsets
    for kk = 1:length(ppg_onsets)
        [~,ppg_find] = min(abs(t_vox-ppg_onsets(kk)));
        temp_response_matrix1(:,:,:,kk) = d1_up(:,:,ppg_find-floor(epoch_pre*srate_epochs):ppg_find+floor(epoch_post*srate_epochs));
        temp_response_matrix2(:,:,:,kk) = d2_up(:,:,ppg_find-floor(epoch_pre*srate_epochs):ppg_find+floor(epoch_post*srate_epochs));
    end
    
    clear d_up d
    
    % output S0:
    response_matrix1(:,:,s,:) = nanmean(temp_response_matrix1,4);
    response_matrix_std1(:,:,s,:) = nanstd(temp_response_matrix1,[],4);
    response_matrix_odd1(:,:,s,:) = nanmean(temp_response_matrix1(:,:,:,1:2:end),4);
    response_matrix_even1(:,:,s,:) = nanmean(temp_response_matrix1(:,:,:,2:2:end),4);
    
    % output T2s:
    response_matrix2(:,:,s,:) = nanmean(temp_response_matrix2,4);
    response_matrix_std2(:,:,s,:) = nanstd(temp_response_matrix2,[],4);
    response_matrix_odd2(:,:,s,:) = nanmean(temp_response_matrix2(:,:,:,1:2:end),4);
    response_matrix_even2(:,:,s,:) = nanmean(temp_response_matrix2(:,:,:,2:2:end),4);
    

    clear temp_response_matrix
end

