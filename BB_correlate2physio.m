function [out_r_map]=BB_correlate2physio(ni,slices)
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

% get the stuff we need:
timing      = BB_get(ni,'timing');
ppg_onsets  = BB_get(ni,'ppg_peaks');
mux_f       = BB_get(ni,'nslices');
srate       = 1/BB_get(ni,'tr');
physio      = BB_get(ni,'physio');

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
    
    out_r_map=zeros(size(ni.data,1),size(ni.data,2),length(slices));

    %%%%% calculate explained variance map

    % cut response_matrix to be one response:
    oneresponse_matrix=temp_resp_mat(:,:,:,t>=0 & t<=1);

    % sample even heartbeats to create a standardized response
    sample_hb=[2:2:size(oneresponse_matrix,3)];
    standard_resp=squeeze(nanmean(oneresponse_matrix(:,:,sample_hb,:),3));
    % figure,plot(squeeze(standard_resp(28,:,:))')

    % predict odd heartbeats
    sample_hb=[1:2:size(oneresponse_matrix,3)];

    % map with predictions
    r_map=NaN(size(oneresponse_matrix,1),size(oneresponse_matrix,2),length(sample_hb));

    for m=1:size(oneresponse_matrix,1)
        disp(['correlation row ' int2str(m) ' of ' int2str(size(oneresponse_matrix,1))])
        for n=1:size(oneresponse_matrix,2)
            for k=1:length(sample_hb) % heartbeats to predict
                % datapoints to predict
                y=squeeze(oneresponse_matrix(m,n,sample_hb(k),:));

                % predictor
                x=squeeze(standard_resp(m,n,:));
                
                % regress
                if ~isempty(y(~isnan(y)))
                    [r,p]=corr(x(~isnan(y)),y(~isnan(y)));
                    r_map(m,n,k)=r;
                end
                clear x y
            end

        end
    end

    out_r_map(:,:,s)=nanmean(r_map,3);
    clear r_map
    

end
