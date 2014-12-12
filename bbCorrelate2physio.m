function [out_r_map]=bbCorrelate2physio(ni,slices)
% function to get correlation with PPG signal
%
% [out_r_map] = BB_correlate2physio(ni)
%
% required input:
% ni: a nifti structure loaded by niftiRead
%
% optional input:
% slices: slices for which to calculate response function, defaults = do
%           all slices
% 
% output:
% out_r_map: voxels X voxels X slices maps with 'correlation' to hearbeat
%
% Wandell Copyright Vistasoft Team, 2013
% Written by Dora 2014

if ~exist('slices','var') % do whole brain
    slices=[1:size(ni.data,3)];
end

% get the nifti stuff we need:
timing      = bbGet(ni,'timing');
mux_f       = bbGet(ni,'super slices');
srate       = 1/bbGet(ni,'tr');

% get physio stuff we need:
physio     = physioCreate('nifti',ni);
ppg_onsets = physioGet(physio,'ppg peaks');

% set epoch times
epoch_pre=0;%sec pre-onset
epoch_post=1;%sec post-onset
step_size=1/srate/mux_f;% in s
srate_epochs=1/step_size;
t=[-epoch_pre:step_size:epoch_post];%s timing for 1 epoch
t_vox=min(timing(:)):step_size:max(timing(:)); % maximal timing accuracy for this scan session

% only use those ppg onsets that have an entire trial before and after
ppg_onsets=ppg_onsets((ppg_onsets-epoch_pre)>0); % get rid of early ones
ppg_onsets=ppg_onsets((ppg_onsets+epoch_post)<max(timing(:))); % get rid of late ones

% create output matrix:
out_r_map=zeros(size(ni.data,1),size(ni.data,2),length(slices));

for s=1:length(slices)
    disp(['slice ' int2str(s) ' of ' int2str(length(slices))])
    sli=slices(s);
    d=squeeze(ni.data(:,:,sli,1:end));

    % get the timing for this slice
    t_sli=timing(sli,:);

    %%%%% get the responses for each voxel, each heartbeat

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
    
    %%%%% calculate explained variance map

    % sample even heartbeats to create a standardized response
    sample_hb=[2:2:size(temp_resp_mat,3)];
    standard_resp=squeeze(nanmean(temp_resp_mat(:,:,sample_hb,:),3));

    % predict odd heartbeats
    sample_hb=[1:2:size(temp_resp_mat,3)];
    
    % map with predictions
    r_map=NaN(size(temp_resp_mat,1),size(temp_resp_mat,2),length(sample_hb));

    for m=1:size(temp_resp_mat,1)
        disp(['correlation row ' int2str(m) ' of ' int2str(size(temp_resp_mat,1))])
        for n=1:size(temp_resp_mat,2)
            for k=1:length(sample_hb) % heartbeats to predict
                % datapoints to predict, this has NaN for the not-sampled
                % points for the trial
                y = squeeze(temp_resp_mat(m,n,sample_hb(k),:));

                % predictor
                x = squeeze(standard_resp(m,n,:));
                
                % regress
                if ~isempty(y(~isnan(y)))
                    r_map(m,n,k)=corr(x(~isnan(y)),y(~isnan(y)));
                end
                clear x y
            end

        end
    end

    out_r_map(:,:,s)=nanmean(r_map,3);
    clear r_map
    

end
