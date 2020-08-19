clear all
close all

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
dDir = '/biac4/wandell/data/BrainBeat/data';

% chdir(dDir)

%% The T2* data are here.  

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

s_nr = 3;
s_info = bb_subs(s_nr);
subj=s_info.subj;


%% Get the anatomicals:

niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

%% Get the MRVenogram:

niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii.gz']));

%%
%% coregister the venogram to the T1:
%%

% if already done, load rotation matrix:
xf_veno=load(fullfile(dDir,subj,veno,[veno_code 'AcpcXform.mat']));

%%
%% now do an SVD on the responses:
%%

% load functional data
scan_nr = 3;
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};
fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
ni = niftiRead(fmri);
% load coregistration matrix:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform.mat']))

% Load the correlation with heartbeat (made with bbCorrelate2physio):
ppgRname = fullfile(dDir,subj,scan,[scanName '_corrPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

% Load the timeseries around PPG peak (made with bbResponse2physio):
ppgTSname = fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT.mat']),'t');

a = reshape(ppgTS.data,[numel(ppgTS.data(:,:,:,1)) length(t)]);
a = a(:,t>=0 & t<=2);
t_sel = t(t>=0 & t<=2);
% a = a-repmat(mean(a,2),1,size(a,2)); % subtract the mean
[u,s,v]=svd(a','econ');
s=diag(s);

% get to cumulative explained variance:
var_explained=cumsum(s.^2) / sum(s.^2)

%%
%% plot a number of components:
%%
nrc_plot=7;

sl_plotx = 26;
sl_plotz = 20;

figure('Position',[0 0 1200 500])
for k=1:nrc_plot
    subplot(3,nrc_plot,k),hold on
%     plot(t_sel,ppgCurve/max(ppgCurve),'Color',[.5 .5 .5])
    plot(t_sel,u(:,k),'b')
    xlim([min(t_sel) max(t_sel)])
    title(['c' int2str(k) ' cumvar ' num2str(var_explained(k),2)])

    whole_brain_v=reshape(v(:,k),[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);

    subplot(3,nrc_plot,nrc_plot+k)
    imagesc(squeeze(whole_brain_v(:,:,sl_plotz)),[-.03 .03])
    axis image
    
    subplot(3,nrc_plot,2*nrc_plot+k)
    im_plot=imrotate(squeeze(whole_brain_v(sl_plotx,:,:)),90);
    imagesc(im_plot,[-.02 .02])
    clear im_plot
    axis image
end

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/2014_12_presentation1/' subj '_' scan '_examples_SVDcurves'])

%%
%% now play around with derivatives
%%

if s_nr==3
    x = 38;
    y = 28;
    z = 20;
end

disp(['r at voxel is ' num2str(ppgR.data(x,y,z))]);
example_curve = double(squeeze(ppgTS.data(x,y,z,:)));

% low pass filter?
srate = 20;
band = 2;
Rp   = 1; Rs = 20; % third order Butterworth
high_p =  band(1)*2/srate;
delta = 0.001*2/srate;
high_s = min(1-delta,high_p+0.1);
[n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
[bf_b,bf_a] = butter(n_band,wn_band,'low');
band_sig    = filtfilt(bf_b,bf_a,example_curve);

figure,
subplot(3,1,1),hold on
plot(t,example_curve)
plot(t,band_sig,'r')

subplot(3,1,2),hold on
plot(t,zeros(size(t)),'k')
plot(t(1:end-1),diff(band_sig),'r')

subplot(3,1,3),hold on
plot(t,zeros(size(t)),'k')
plot(t(1:end-2),diff(diff(band_sig)),'r')

%%
%% take a slice
%%

a = double(reshape(ppgTS.data(:,:,:,:),[numel(ppgTS.data(:,:,:,1)) length(t)]));

% low pass filter?
srate = 20;
band = 2;
Rp   = 1; Rs = 20; % third order Butterworth
high_p =  band(1)*2/srate;
delta = 0.001*2/srate;
high_s = min(1-delta,high_p+0.1);
[n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
[bf_b,bf_a] = butter(n_band,wn_band,'low');
for k=1:size(a,1)
    a(k,:) = filtfilt(bf_b,bf_a,a(k,:));
end
a = a(:,t>=0 & t<=2);
a = diff(a,2,2);
t_sel = t(t>=0 & t<=2);
t_sel = t_sel(1:end-2);
% a = a-repmat(mean(a,2),1,size(a,2)); % subtract the mean
[u,s,v]=svd(a','econ');
s=diag(s);

% get to cumulative explained variance:
var_explained=cumsum(s.^2) / sum(s.^2)

%%
nrc_plot=7;

sl_plotx = 26;
sl_plotz = 20;

figure('Position',[0 0 1200 500])
for k=1:nrc_plot
    subplot(3,nrc_plot,k),hold on
%     plot(t_sel,ppgCurve/max(ppgCurve),'Color',[.5 .5 .5])
    plot(t_sel,u(:,k),'b')
    xlim([min(t_sel) max(t_sel)])
    title(['c' int2str(k) ' cumvar ' num2str(var_explained(k),2)])

    whole_brain_v=reshape(v(:,k),[size(ppgTS.data,1) size(ppgTS.data,2) size(ppgTS.data,3)]);

    subplot(3,nrc_plot,nrc_plot+k)
    imagesc(squeeze(whole_brain_v(:,:,sl_plotz)),[-.03 .03])
    axis image
    
    subplot(3,nrc_plot,2*nrc_plot+k)
    im_plot=imrotate(squeeze(whole_brain_v(sl_plotx,:,:)),90);
    imagesc(im_plot,[-.02 .02])
    clear im_plot
    axis image
end

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/2014_12_presentation1/' subj '_' scan '_examples_SVDcurves'])

