function [response_matrix,t,response_matrix_odd,response_matrix_even,response_matrix_std,response_matrix_sterr] = ...
    bbResponse2physio(ni,slices,varargin)
% function to get brain response after the peak of a heartbeat (or
% respiration later)
%
% [response_matrix,t,response_matrix_odd,response_matrix_even] = bbResponse2physio(ni,slices)
% 
% ni: a nifti structure loaded by niftiRead
%
% optional inputs:
% slices: slices for which to calculate response function, if no input, do
% all slices
%
% [response_matrix,t,response_matrix_odd,response_matrix_even] = bbResponse2physio(ni,slices)

% DHermes, Copyright Vistasoft Team 2014

if ~exist('slices','var') % do whole brain
    slices = [1:size(ni.data,3)];
elseif exist('slices','var') && isempty(slices)% do whole brain
    slices = [1:size(ni.data,3)];
end

% get the nifti stuff we need:
timing = bbGet(ni,'timing'); % timing per slice
mux_f = bbGet(ni,'super slices');
srate = 1/bbGet(ni,'tr');

% get physio stuff we need:
physio = physioCreate('nifti',ni);
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
step_size = 1/srate/mux_f; % in s
srate_epochs = 1/step_size;
t = [-epoch_pre:step_size:epoch_post]; % s timing for 1 epoch

% maximal timing accuracy for this scan session
t_vox = min(timing(:)):step_size:(max(timing(:)+step_size)); 

% only use those ppg onsets that have an entire trial before and after
ppg_onsets = ppg_onsets((ppg_onsets-epoch_pre)>1); % get rid of early ones, and the first scans
ppg_onsets = ppg_onsets((ppg_onsets+epoch_post)<max(timing(:))); % get rid of late ones

% initiate output matrix to fill:
response_matrix = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));
response_matrix_std = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));
response_matrix_sterr = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));
response_matrix_odd = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));
response_matrix_even = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));

for ss = 1:length(slices)
    disp(['slice ' int2str(ss) ' of ' int2str(length(slices))])
    sli = slices(ss);
    d = squeeze(ni.data(:,:,sli,1:end));

    % demean and express in percent modulation 
    d_norm = reshape(d,[size(d,1) * size(d,2), size(d,3)]);
    points_use = 4:size(d_norm,2); % do not use the first three scans
    for kk = 1:size(d_norm,1) % run through voxels in this slice
        % do not use first scans:
        x = points_use;
        y = d_norm(kk,points_use);
        
        % detrend
        p = polyfit(x,y,1);    

        % percent modulation
        mean_factor = mean(d_norm(kk,points_use)); % mean
        d_norm(kk,:) = (d_norm(kk,:) - (p(1)*[1:size(d_norm,2)] + p(2)))./mean_factor;
    end
    d = reshape(d_norm,[size(d,1), size(d,2), size(d,3)]);
    clear d_norm x y

    % create a set of NaNs at the upsampled time rate 
    d_up = single(NaN(size(d,1),size(d,2),length(t_vox))); % initiate with NaNs

    % get the timing for this slice (times where we have data)
    t_sli = timing(sli,:);

    % fill the data into these times:
    % get the indices in upsampled time at which data were collected
    d_up(:,:,ismember(round(t_vox,3),round(t_sli,3))) = d;

    % temporary response matrix for 1 slice: voxel X voxel X time X PPGonset
    temp_response_matrix = single(zeros(size(ni.data,1),size(ni.data,2),length(t),length(ppg_onsets)));
    % run through all ppg onsets
    for k = 1:length(ppg_onsets)
        [~,ppg_find] = min(abs(t_vox-ppg_onsets(k)));
        temp_response_matrix(:,:,:,k) = d_up(:,:,ppg_find-floor(epoch_pre*srate_epochs):ppg_find+floor(epoch_post*srate_epochs));
    end
    
    clear d_up d
    
    % output:
    response_matrix(:,:,ss,:) = mean(temp_response_matrix,4,'omitnan');
    response_matrix_std(:,:,ss,:) = nanstd(temp_response_matrix,[],4);
    response_matrix_odd(:,:,ss,:) = nanmean(temp_response_matrix(:,:,:,1:2:end),4);
    response_matrix_even(:,:,ss,:) = nanmean(temp_response_matrix(:,:,:,2:2:end),4);
    
    % create standard error output as well, but first calculate number of
    % NaNs in each timepoint
    nr_nonnans = sum(~isnan(temp_response_matrix(1,1,:,:)),4); % just need for 1 voxel, similar for other voxels in slice
    response_matrix_sterr(:,:,ss,:) = nanstd(temp_response_matrix,[],4)./sqrt(repmat(nr_nonnans,size(temp_response_matrix,1),size(temp_response_matrix,2),1));

    clear temp_response_matrix
end

