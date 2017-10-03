function [response_matrix,t,response_matrix_odd,response_matrix_even] = bbResponse2physio(ni,slices,varargin)
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

% Written by Dora, Copyright Vistasoft Team 2014

if ~exist('slices','var') % do whole brain
    slices=[1:size(ni.data,3)];
elseif exist('slices','var') && isempty(slices)% do whole brain
    slices=[1:size(ni.data,3)];
end

% get the nifti stuff we need:
timing=bbGet(ni,'timing');
mux_f=bbGet(ni,'super slices');
srate=1/bbGet(ni,'tr');

% get physio stuff we need:
physio     = physioCreate('nifti',ni);
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
t_vox = min(timing(:)):step_size:max(timing(:)); % maximal timing accuracy for this scan session

% only use those ppg onsets that have an entire trial before and after
ppg_onsets = ppg_onsets((ppg_onsets-epoch_pre)>1); % get rid of early ones, and the first scans
ppg_onsets = ppg_onsets((ppg_onsets+epoch_post)<max(timing(:))); % get rid of late ones

% initiate output matrix to fill:
response_matrix = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));
response_matrix_odd = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));
response_matrix_even = single(zeros(size(ni.data,1),size(ni.data,2),length(slices),length(t)));

for s = 1:length(slices)
    disp(['slice ' int2str(s) ' of ' int2str(length(slices))])
    sli = slices(s);
    d = squeeze(ni.data(:,:,sli,1:end));

    % demean and express in percent modulation 
    d_norm=reshape(d,[size(d,1) * size(d,2), size(d,3)]);
    points_use=4:size(d_norm,2); % do not use the first couple of scans
    for k=1:size(d_norm,1)
        % do not use first scans:
        x = points_use;
        y = d_norm(k,points_use);
        
        % detrend
        p = polyfit(x,y,1);    

        % percent modulation
        mean_factor = mean(d_norm(k,points_use)); % mean
        d_norm(k,:) = (d_norm(k,:) - (p(1)*[1:size(d_norm,2)] + p(2)))./mean_factor;
    end
    d=reshape(d_norm,[size(d,1), size(d,2), size(d,3)]);
    clear d_norm

    % get the timing for this slice
    t_sli = timing(sli,:);

    % get the upsampled (NaN) of this slice
    d_up=single(NaN(size(d,1),size(d,2),length(t_vox))); % initiate with NaNs
    d_up(:,:,ismember(round(t_vox*mux_f*srate),round(t_sli*mux_f*srate)))=d;

    temp_response_matrix = single(zeros(size(ni.data,1),size(ni.data,2),length(t),length(ppg_onsets)));
    % run through all ppg onsets
    for k=1:length(ppg_onsets); 
        [~,ppg_find]=min(abs(t_vox-ppg_onsets(k)));
        temp_response_matrix(:,:,:,k)=d_up(:,:,ppg_find-round(epoch_pre*srate_epochs):ppg_find+round(epoch_post*srate_epochs));
    end
    
    clear d_up d
    
    % output:
    response_matrix(:,:,s,:) = nanmean(temp_response_matrix,4);
    response_matrix_odd(:,:,s,:) = nanmean(temp_response_matrix(:,:,:,1:2:end),4);
    response_matrix_even(:,:,s,:) = nanmean(temp_response_matrix(:,:,:,2:2:end),4);
    clear temp_response_matrix
end

