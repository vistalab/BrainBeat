function [response_matrixTE1,t,response_matrix_oddTE1,response_matrix_evenTE1,response_matrix_stdTE1,...
    response_matrixTE2,response_matrix_oddTE2,response_matrix_evenTE2,response_matrix_stdTE2] = ...
        bbResponse2physioME_te1te2(ni1,ni2,slices,varargin)
%
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
% Returns Log(S0) and RR or T2*?
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
physio = physioCreate('nifti',ni1);
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

% initiate output matrices to fill for TE1 and TE2
response_matrixTE1 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_stdTE1 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_oddTE1 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_evenTE1 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));

response_matrixTE2 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_stdTE2 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_oddTE2 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_evenTE2 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));

for s = 1:length(slices)
    
    disp(['slice ' int2str(s) ' of ' int2str(length(slices))])
    sli = slices(s);
    
    % get the log for this slice, so we can solve the linear form
    d1 = log(squeeze(ni1.data(:,:,sli,1:end)));
    d2 = log(squeeze(ni2.data(:,:,sli,1:end)));
    
    % detrend now
    sigte1_norm = reshape(d1,[size(d1,1) * size(d1,2), size(d1,3)]);
    sigte2_norm = reshape(d2,[size(d2,1) * size(d2,2), size(d2,3)]);
    points_use = 4:size(sigte1_norm,2); % do not use the first couple of scans
    for kk = 1:size(sigte1_norm,1) % voxels in scans
        % do not use first scans:
        x = points_use;
        
        % TE1 & TE2 signals
        y = (sigte1_norm(kk,points_use) + sigte2_norm(kk,points_use))/2;
        % detrend mean of TE1 and TE2
        p = polyfit(x,y,1);    
        % apply to TE1
        sigte1_norm(kk,:) = sigte1_norm(kk,:) - (p(1)*[1:size(sigte1_norm,2)]);
        % apply to TE2
        sigte2_norm(kk,:) = sigte2_norm(kk,:) - (p(1)*[1:size(sigte2_norm,2)]);
        
%         % alternatively, one could perhaps detrend each echo, but means are
%         % subtracted, so S0 would be lost in this case...
%         % apply to TE1
%         sigte1_norm(kk,:) = detrend(sigte1_norm(kk,:));
%         % apply to TE2
%         sigte2_norm(kk,:) = detrend(sigte2_norm(kk,:));
    end
        
    d1 = reshape(sigte1_norm,[size(d1,1), size(d1,2), size(d1,3)]);
    d2 = reshape(sigte2_norm,[size(d2,1), size(d2,2), size(d2,3)]);
    clear sigte1_norm sigte2_norm

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
    temp_response_matrixTE1 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(t),length(ppg_onsets)));
    temp_response_matrixTE2 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(t),length(ppg_onsets)));
    % run through all ppg onsets
    for kk = 1:length(ppg_onsets)
        [~,ppg_find] = min(abs(t_vox-ppg_onsets(kk)));
        temp_response_matrixTE1(:,:,:,kk) = d1_up(:,:,ppg_find-floor(epoch_pre*srate_epochs):ppg_find+floor(epoch_post*srate_epochs));
        temp_response_matrixTE2(:,:,:,kk) = d2_up(:,:,ppg_find-floor(epoch_pre*srate_epochs):ppg_find+floor(epoch_post*srate_epochs));
    end
    
    clear d_up d
        
    % output TE1:
    response_matrixTE1(:,:,s,:) = nanmean(temp_response_matrixTE1,4);
    response_matrix_stdTE1(:,:,s,:) = nanstd(temp_response_matrixTE1,[],4);
    response_matrix_oddTE1(:,:,s,:) = nanmean(temp_response_matrixTE1(:,:,:,1:2:end),4);
    response_matrix_evenTE1(:,:,s,:) = nanmean(temp_response_matrixTE1(:,:,:,2:2:end),4);

    % output TE2:
    response_matrixTE2(:,:,s,:) = nanmean(temp_response_matrixTE2,4);
    response_matrix_stdTE2(:,:,s,:) = nanstd(temp_response_matrixTE2,[],4);
    response_matrix_oddTE2(:,:,s,:) = nanmean(temp_response_matrixTE2(:,:,:,1:2:end),4);
    response_matrix_evenTE2(:,:,s,:) = nanmean(temp_response_matrixTE2(:,:,:,2:2:end),4);

    clear temp_response_matrix
end

