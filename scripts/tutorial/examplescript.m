clear all
close all

data_path='/biac4/wandell/data/BrainBeat/data/';
subj_name='20140911_1411';
scan_name='3_1_mux8fov4_r1_35s_3mm';
fmri_data=fullfile(data_path,subj_name,scan_name,'7944_3_1.nii.gz');
scan_name='8_1_T1w_1mm_sag';
anat_data=fullfile(data_path,subj_name,scan_name,'/7944_8_1.nii.gz');

ni=niftiRead(fmri_data);

ni=niftiRead(anat_data);

%%

% get slice times:
cd /green/
nslices = BB_get(ni,'mux');
 

%%

figure
for k=1:size(ni.data,4)
imagesc(double(ni.data(:,:,30,k))); axis image;
colormap(hot)
pause(.35)

end


%% play around with latency maps

[out_p_map]=bbPhase2physio(ni,20);

%% try SVD on the whole brain curves:

niRespPPG = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse.nii.gz']));
load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT.mat']) ,'t') % load time t

a = reshape(niRespPPG.data,[numel(niRespPPG.data(:,:,:,1)) length(t)]);
a = a-repmat(mean(a,2),1,size(a,2)); % subtract the mean
[u,s,v]=svd(a','econ');
s=diag(s);

% check whether this needs to be squared:
var_explained=cumsum(s.^2) / sum(s.^2)

%%

% pick a slice
sl_plot=28;

nrc_plot=10;
figure
for k=1:nrc_plot
    subplot(2,nrc_plot,k),hold on
%     plot(t,ppgCurve/max(ppgCurve),'Color',[.5 .5 .5])
    plot(t,u(:,k),'b')
    xlim([min(t) max(t)])
    title(['comp ' int2str(k)])

    whole_brain_v=reshape(v(:,k),[size(niRespPPG.data,1) size(niRespPPG.data,2) size(niRespPPG.data,3)]);

    subplot(2,nrc_plot,nrc_plot+k)
    imagesc(squeeze(whole_brain_v(sl_plot,:,:)),[-.03 .03])
end


%% % try SVD on the response curves

a=reshape(standard_resp,[56*56,25]);
% figure,plot(a')
a=a-repmat(mean(a,2),1,size(a,2));

figure,plot(a')

[u,s,v]=svd(a','econ');

% now make a plot of a number of components:
nrc_plot=10;
figure
for k=1:nrc_plot
    subplot(2,nrc_plot,k),hold on
    plot(t,ppgCurve/max(ppgCurve),'Color',[.5 .5 .5])
    plot(t,u(:,k),'b')
    title(['comp ' int2str(k)])
    
    subplot(2,nrc_plot,nrc_plot+k)
    imagesc(reshape(v(:,k),[56,56]),[-.1 .1])
end


