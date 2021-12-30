function [response_matrixS0,t,response_matrix_oddS0,response_matrix_evenS0,response_matrix_stdS0,...
    response_matrixT2s,response_matrix_oddT2s,response_matrix_evenT2s,response_matrix_stdT2s] = ...
        bbResponse2physioME_s0t2s(ni1,ni2,slices,varargin)
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

% initiate output matrices to fill for TE1, TE2, S0 and T2*
response_matrixS0 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_stdS0 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_oddS0 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_evenS0 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));

response_matrixT2s = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_stdT2s = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_oddT2s = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));
response_matrix_evenT2s = single(zeros(size(ni1.data,1),size(ni1.data,2),length(slices),length(t)));

for s = 1:length(slices)
    
    disp(['slice ' int2str(s) ' of ' int2str(length(slices))])
    sli = slices(s);
    
    % get the log for this slice, so we can solve the linear form
    d1 = log(squeeze(ni1.data(:,:,sli,1:end)));
    d2 = log(squeeze(ni2.data(:,:,sli,1:end)));

    % get echo times in units of seconds
    te1 = 0.001*bbGet(ni1,'te');
    te2 = 0.001*bbGet(ni2,'te');
    
    % calculate T2*
    rr = (d2 - d1)/(te2 - te1);
    t2s = -1./rr;
    
    % calculate log S0, exponentiate later?
    lns0 = d1 - (d2 - d1)*(te1/(te2 - te1));
    
%     % plot to check for subject 4, multi echo run 1
%     if s==20
%         fn_MakeSuppFigME(d1,d2,lns0,rr,t2s,srate,ppg_onsets)
%     end

    % detrend now
    lns0_norm = reshape(lns0,[size(lns0,1) * size(lns0,2), size(lns0,3)]);
    t2s_norm = reshape(t2s,[size(t2s,1) * size(t2s,2), size(t2s,3)]);
    points_use = 4:size(lns0_norm,2); % do not use the first couple of scans
    for kk = 1:size(lns0_norm,1)
        % do not use first scans:
        x = points_use;
        
        % S0
        y = lns0_norm(kk,points_use);       
        % fit line
        p = polyfit(x,y,1);           
        % linear detrend only, keep mean
        lns0_norm(kk,:) = lns0_norm(kk,:) - (p(1)*[1:size(lns0_norm,2)]);

        % T2*
        y = t2s_norm(kk,points_use);       
        % detrend
        p = polyfit(x,y,1);    
        % linear detrend only, keep mean
        t2s_norm(kk,:) = t2s_norm(kk,:) - (p(1)*[1:size(t2s_norm,2)]);
        
    end
    
%     % we could exponentiate to get S0
%     d1_norm = exp(d1_norm);
    
    d1 = reshape(lns0_norm,[size(d1,1), size(d1,2), size(d1,3)]);
    d2 = reshape(t2s_norm,[size(d2,1), size(d2,2), size(d2,3)]);
    clear lns0_norm t2s_norm

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
    temp_response_matrixS0 = single(zeros(size(ni1.data,1),size(ni1.data,2),length(t),length(ppg_onsets)));
    temp_response_matrixT2s = single(zeros(size(ni1.data,1),size(ni1.data,2),length(t),length(ppg_onsets)));
    % run through all ppg onsets
    for kk = 1:length(ppg_onsets)
        [~,ppg_find] = min(abs(t_vox-ppg_onsets(kk)));
        temp_response_matrixS0(:,:,:,kk) = d1_up(:,:,ppg_find-floor(epoch_pre*srate_epochs):ppg_find+floor(epoch_post*srate_epochs));
        temp_response_matrixT2s(:,:,:,kk) = d2_up(:,:,ppg_find-floor(epoch_pre*srate_epochs):ppg_find+floor(epoch_post*srate_epochs));
    end
    
    clear d_up d
    
    % output S0:
    response_matrixS0(:,:,s,:) = nanmean(temp_response_matrixS0,4);
    response_matrix_stdS0(:,:,s,:) = nanstd(temp_response_matrixS0,[],4);
    response_matrix_oddS0(:,:,s,:) = nanmean(temp_response_matrixS0(:,:,:,1:2:end),4);
    response_matrix_evenS0(:,:,s,:) = nanmean(temp_response_matrixS0(:,:,:,2:2:end),4);
    
    % output T2s:
    response_matrixT2s(:,:,s,:) = nanmean(temp_response_matrixT2s,4);
    response_matrix_stdT2s(:,:,s,:) = nanstd(temp_response_matrixT2s,[],4);
    response_matrix_oddT2s(:,:,s,:) = nanmean(temp_response_matrixT2s(:,:,:,1:2:end),4);
    response_matrix_evenT2s(:,:,s,:) = nanmean(temp_response_matrixT2s(:,:,:,2:2:end),4);

    clear temp_response_matrix
end

