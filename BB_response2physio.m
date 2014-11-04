function [response_matrix,t] = BB_response2physio(ni,slices)
% function to get brain response after the peak of a heartbeat (or
% respiration later)
%
% [response_matrix,t] = BB_response2physio(ni)
% 
% ni: a nifti structure loaded by niftiRead
%
% optional inputs:
% slices: slices for which to calculate response function, if no input, do
% all slices
% 
% 
% Written by Dora, Copyright Vistasoft Team 2014

if ~exist('slices','var') % do whole brain
    slices=[1:size(ni.data,3)];
end

% get the stuff we need:
timing=bbGet(ni,'timing');
ppg_onsets=bbGet(ni,'ppg_peaks');
mux_f=bbGet(ni,'super slices');
srate=1/bbGet(ni,'tr');

% set epoch times
epoch_pre=.5;%sec pre-onset
epoch_post=2;%sec post-onset
step_size=1/srate/mux_f;% in s
srate_epochs=1/step_size;
t=[-epoch_pre:step_size:epoch_post];%s timing for 1 epoch
t_vox=min(timing(:)):step_size:max(timing(:)); % maximal timing accuracy for this scan session

% initiate output matrix to fill:
response_matrix=zeros(size(ni.data,1),size(ni.data,2),length(t),length(slices));

% only use those ppg onsets that have an entire trial before and after
ppg_onsets=ppg_onsets((ppg_onsets-epoch_pre)>0); % get rid of early ones
ppg_onsets=ppg_onsets((ppg_onsets+epoch_post)<max(timing(:))); % get rid of late ones

for s=1:length(slices)
    disp(['slice ' int2str(s) ' of ' int2str(length(slices))])
    sli=slices(s);
    d=squeeze(ni.data(:,:,sli,1:end));

    % get the timing for this slice
    t_sli=timing(sli,:);

    % walk through voxels in the slice:
    temp_resp_mat=NaN(size(d,1),size(d,2),length(ppg_onsets),length(t));
    for m=1:size(d,1)
        for n=1:size(d,2)
            % get the upsampled (NaN) of this voxel
            s_vox=squeeze(d(m,n,:)); % signal
            s_vox_up=NaN(length(t_vox),1); % initiate with NaNs
            s_vox_up(ismember(round(t_vox*mux_f*srate),round(t_sli*mux_f*srate)))=s_vox;

            for k=1:length(ppg_onsets); % run through all ppg onsets
                [~,ppg_find]=min(abs(t_vox-ppg_onsets(k)));
                temp_resp_mat(m,n,k,:)=s_vox_up(ppg_find-round(epoch_pre*srate_epochs):ppg_find+round(epoch_post*srate_epochs));
            end
            clear s_vox
        end
    end
    
    response_matrix(:,:,:,s)=squeeze(nanmean(temp_resp_mat,3));
    clear temp_resp_mat
end