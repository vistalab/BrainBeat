function [out_p_map]=bbPhase2physio(ni,slices)
% function to get correlation with PPG signal
%
% [out_p_map] = BB_phase2physio(ni)
%
% required input:
% ni: a nifti structure loaded by niftiRead
%
% optional input:
% slices: slices for which to calculate response function, defaults = do
%           all slices
% 
% output:
% out_p_map: voxels X voxels X slices maps with 'phase' to hearbeat
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
ppgSrate = physioGet(physio,'ppg srate');
ppgCurve = physioGet(physio,'ppg curve'); % get one PPG curve
% resample at brain data epoch srate
ppgCurve=resample(ppgCurve,round(1/max(diff(sort(timing(:))))),ppgSrate);
ppgRate = physioGet(physio,'ppg rate'); %heart rate

% set epoch times
epoch_pre=0;%sec pre-onset
epoch_post=1/ppgRate;%sec post-onset
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
    % hearbeat_response averaged
    standard_resp=squeeze(nanmean(temp_resp_mat,3));
    

    %%%% calculate 'phase' map
    % create output:
    out_p_map=zeros(size(ni.data,1),size(ni.data,2),length(slices));

    % map with phase (peak latency) predictions
    p_map=NaN(size(standard_resp,1),size(standard_resp,2));

    for m=1:size(standard_resp,1)
        disp(['phase row ' int2str(m) ' of ' int2str(size(standard_resp,1))])
        for n=1:size(standard_resp,2)
            % response for this voxel
            x=squeeze(standard_resp(m,n,:));
            
            % compute the cross-correlation
            [x_corr,x_lags] = crosscorr(x,ppgCurve,length(x)-1);
            % find the maximum (absolute)
            [~,max_x_corr_i] = max(abs(x_corr));
            % find the lag in time
            t_lag = x_lags(max_x_corr_i) * step_size;
            
            % COMMENT: I still do not really like this latency calculation,
            % does not seem to capture a lot of information
            % TODO: find out whether it is positive or negative
            
            p_map(m,n) = t_lag;
            
            clear x x_f x_an phase_shift x_peaks
        end
    end

    out_p_map(:,:,s)=p_map;
    clear p_map
end