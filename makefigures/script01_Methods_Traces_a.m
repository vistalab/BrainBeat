clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%
%

%% Base data directory on a Mac mounting biac4 (wandell's machine)

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';
% chdir(dDir)

%% The T2* data are here.  

% Select a subject and scan nummer
s_nr = 2;
scan_nr = 1;

subs = bb_subs(s_nr);
subj = subs.subj;
scan = subs.scan{scan_nr};
scanName = subs.scanName{scan_nr};
fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
if ~exist(fmri,'file')
    clear ni
    error('filename %s does not exist',fmri)
end
ni = niftiRead(fmri);

% load coregistration matrix for the functionals:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

%% Anatomicals

anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii']);
niAnatomy = niftiRead(anat);

%% quick overlay between functionals and anatomy
sliceThisDim = 1;
if s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [1,10,-20];
%     curPos = [-29,-14,-50]; % 03
%     curPos = [-3,30,-43]; % 03
    curPos = [1,10,-20];
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [0,4,38];
%     curPos = [0,4,38];
end
niFunc = ni;

% niFunc.data = ni.data(:,:,:,1); % overlay the first functional - more structure visible
% bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
% title('First functional on anatomy')
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/checkCoreg/' subj '_' scan '_Func1onAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

niFunc.data = mean(ni.data(:,:,:,5:end),4); % overlay the mean functional
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
title('Mean functional on anatomy')
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/checkCoreg/' subj '_' scan '_MeanFuncOnAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])


%% Plot physiology data

figure('Position',[0 0 600 600])

% get data
physio      = physioCreate('nifti',ni);
ppgData     = physioGet(physio,'ppg data');
respData    = physioGet(physio,'resp data');

% plot PPG peaks on PPG signal
subplot(3,4,1:3)
ppgPeaks    = physioGet(physio,'ppg peaks');
tPPG = physioGet(physio,'ppg sample times');
srate = physioGet(physio,'ppg srate');
ppgPeakSamples = round(ppgPeaks*srate);
plot(ppgPeaks,ppgData(ppgPeakSamples),'ro',tPPG,ppgData,'k-');
xlabel('time (s)'); grid on
title('PPG')
xlim([0 10])

% plot one PPG curve
subplot(3,4,4)
ppgCurve = physioGet(physio,'ppg curve');
ppgSrate = physioGet(physio,['PPGsrate']); %  = physio.ppg.srate
plot((1:length(ppgCurve))/ppgSrate,ppgCurve,'k','LineWidth',2);
xlim([0 max((1:length(ppgCurve))/ppgSrate)])
xlabel('time (s)')
% get PPG rate
ppgRate = physioGet(physio,'ppg rate');
title(['PPG rate = ' num2str(ppgRate)])

% plot RESP peaks on RESP signal
subplot(3,4,5:7)
respPeaks = physioGet(physio,'resp peaks');
tResp = physioGet(physio,'resp sample times');
srate = physioGet(physio,'resp srate');
respPeakSamples = round(respPeaks*srate);
plot(respPeaks,respData(respPeakSamples),'ro',tResp,respData,'k-');
xlabel('time (s)'); grid on
xlim([0 50])
title('RESP')

% plot one RESP curve
subplot(3,4,8)
respCurve = physioGet(physio,'resp curve');
respSrate = physioGet(physio,['RESPsrate']);
plot((1:length(respCurve))/respSrate,respCurve,'k','LineWidth',2);
xlabel('time (s)')
xlim([0 max((1:length(respCurve))/respSrate)])
% get RESP rate
respRate = physioGet(physio,'resp rate');
title(['RESP rate = ' num2str(respRate)])

% plot PPG signal with scan onsets
subplot(3,1,3),hold on
plot(ppgPeaks,ppgData(ppgPeakSamples),'ro',tPPG,ppgData,'k-');
plot(tPPG(physio.ppg.scan_onset==1),0,'b*')
title('PPG signal (black) and scan onsets (blue)')
xlim([0 10])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-depsc',[dDir './figures/physio/' subj '_' scan '_physioTrace'])
print('-painters','-r300','-dpng',[dDir './figures/physio/' subj '_' scan '_physioTrace'])

%% Plot MRI data

in_data = 'PPG';

sliceThisDim = 3;
if s_nr == 2
    imDims = [-90 -120 -120; 90 130 90];
%     curPos = [1,10,-20];
%     curPos = [-29,-14,-50]; % 03
%     curPos = [-3,30,-43]; % 03
    curPos = [1,10,-30];
elseif s_nr == 3
    imDims = [-90 -120 -100; 90 130 110];
    curPos = [1,4,38];
end

% load time series and associated time
ppgTSname = fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponse.nii']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponseT.mat']),'t');

% load correlation between even and odd scans for colors of timeseries
ppgRname = fullfile(dDir,subj,scan,[scanName '_corr' in_data '.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

%%%% Overlay 1: functionals and anatomy
niFunc = ni;
niFunc.data = mean(ni.data(:,:,:,5:end),4); % overlay the mean functional
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)

%%%% Overlay 2: timeseries and anatomy
ppgTSplot = ppgTS;
if isequal(in_data,'PPG')
    ppgTSplot.data(:,:,:,t<-.1 | t>1)=[]; % plot these times from curve
elseif isequal(in_data,'RESP')
    ppgTSplot.data(:,:,:,t<-.2 | t>4)=[]; % plot these times from curve
end
niColor = ppgR; % use R2 map to scale later: r2_scale=sqrt(imgSlice2.^2);
niColor.data = niColor.data.^2; %r^2 - squared correlation coefficient between even and odd beats
bbOverlayTimeseriesAnat(ppgTSplot,niColor,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)

clear niColor ppgTSplot

%% get MRI responses after PPG peak

sl_plot = 20;
[response_matrix,t,temp_resp_mat] = bbResponse2physio(ni,sl_plot);

% pick one voxel that we know is good:
m=28; n=35;
a=squeeze(temp_resp_mat(m,n,:,:));

figure('Position',[0 0 600 400])
subplot(1,2,1),hold on
plot(t,a(1,:),'r.')
plot(t,a(2,:),'g.')
plot(t,a(3,:),'b.')
xlim([min(t) max(t)])
title('heartbeat 1, 2 and 3')

subplot(1,2,2),hold on
plot(t,a','k.')
plot(t,nanmean(a(1:3:end,:),1),'b','LineWidth',2)
plot(t,nanmean(a,1),'r','LineWidth',2)
xlim([min(t) max(t)])
title('red - all; blue - 1/3 of data')

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/2014_12_presentation1/' subj '_' scan '_brain_ppg1'])



%%
%%
%%
%%
%% OLD
%%
%%
%%
%%
 

%% ONLY MOVIE of 1 slice

% now run for the whole brain:
[response_matrix,t] = bbResponse2physio(ni);

% safe the response matrix as mat
save([data_path subj_name '/' scan_name '/' b '_PPGtrigResponse'],'response_matrix','t')

% safe the response matrix as nifti and save  time t
ni1=ni;
ni1.data=response_matrix;
[~,b]=fileparts(ni.fname);
[~,b]=fileparts(b);
ni1.fname=[b '_PPGtrigResponse'];
niftiWrite(ni1,[data_path subj_name '/' scan_name '/' ni1.fname])
clear ni1
save([data_path subj_name '/' scan_name '/' b '_PPGtrigResponseT'],'t')


%% quick figure for one slice, one row
figure('Position',[0 0 500 800])
plot(t,squeeze(response_matrix(28,:,:,20))) % 28/35

% get ppg response curve 
phys_srate=physio.ppg.srate;
signal=physio.ppg.data;
epoch_pre=-t(1)*phys_srate;
epoch_post=t(end)*phys_srate;
ppg_t=(-round(epoch_pre)+1:epoch_post)/phys_srate;
phys_epochs1=zeros(length(ppg_onsets)-6,length(ppg_t));
for k=2:length(ppg_onsets)-5
    sample_ppg=round(ppg_onsets(k)*phys_srate);
    phys_epochs1(k,:)=signal(sample_ppg-round(epoch_pre)+1:sample_ppg+(epoch_post));
end
% get average ppg response
ppg_resp=mean(phys_epochs1,1);

%% make the movie:

% pick a slice:
sl_plot=28;
a=squeeze(response_matrix(sl_plot,:,:,:)); % third dimension is time

%%

nFrames=length(t);
mov(1:nFrames)=struct('cdata',[],'colormap',[]);

figure('Position',[0 0 400 400])
subplot(2,1,2),hold on
plot(ppg_t,ppg_resp-mean(ppg_resp));

a_meansub=a-repmat(mean(a,3),[1 1 size(a,3)]);

for k=1:size(a,3)
    subplot(2,2,1)
    imagesc(flipud(a(:,:,k)'),[0 max(a(:))])
%     imagesc(a(:,:,k),[0 max(a(:))])
    axis image
    title(['t=' num2str(t(k))])
    
    subplot(2,2,2)
    imagesc(flipud(a_meansub(:,:,k)'),[-2 2])
%     imagesc(a_meansub(:,:,k),[-2 2])
    axis image
    title(['t=' num2str(t(k))])

    subplot(2,1,2),hold on
    xlim([t(1) t(end)])
    plot(t(k),0,'r.')
    
    pause(1)
    mov(k)=getframe(gcf);
end
close gcf

% [~,b]=fileparts(ni.fname);
% [~,b]=fileparts(b);
% movie2avi(mov,fullfile(data_path,subj_name,scan_name,[b '_PPGmov_slice' int2str(sl_plot) '.avi']),'compression','None','fps',1)

%% make figure with curves
sl_plot=[28 1];
% sl_plot=[28 2];
% sl_plot=[20 3];

if sl_plot(2)==1
    plotting_scan1=squeeze(ni.data(sl_plot(1),:,:,1))';
    a=squeeze(response_matrix(sl_plot(1),:,:,t>-.2 & t<1)); % fourth dimension is time
    r2_scale=sqrt(squeeze(out_r_map(sl_plot(1),:,:)).^2);
elseif sl_plot(2)==2
    plotting_scan1=squeeze(ni.data(:,sl_plot(1),:,1))';
    a=squeeze(response_matrix(:,sl_plot(1),:,t>-.2 & t<1)); % fourth dimension is time
    r2_scale=sqrt(squeeze(out_r_map(:,sl_plot(1),:)).^2);
elseif sl_plot(2)==3
    plotting_scan1=squeeze(ni.data(:,:,sl_plot(1),1))';
    a=squeeze(response_matrix(:,:,sl_plot(1),t>-.2 & t<1)); % fourth dimension is time
    r2_scale=sqrt(squeeze(out_r_map(:,:,sl_plot(1))).^2);
end

f=figure;
cm=colormap(jet);
close(f)

figure('Position',[0 0 1200 800]),hold on
imagesc(plotting_scan1)
colormap gray

for k=1:size(a,1)
    for m=1:size(a,2)
        
        % subtract the mean:
        r_curve=squeeze(a(k,m,:)-mean(a(k,m,:)));
        
        % scale by strength of correlation with heartbeat:
%         r_curve=r_curve/max(abs(r_curve(:)));
%         r_curve=r_curve*r2_scale(k,m);
%         plot(k+[0:size(a,3)-1]/size(a,3)-.5,m+r_curve,'r')

        % just scale by the maximum and color with correlation strength
        if max(r_curve)>.8
            r_curve=.8*r_curve/max(r_curve);
        end
        c_use=cm(1+floor(r2_scale(k,m)*64),:);
        plot(k+[0:size(a,3)-1]/size(a,3)-.5,m+r_curve,'Color',c_use)
        
        plot(k-.5,m,'.','Color',[.5 .5 .5])
    end
end
axis image

set(gcf,'PaperPositionMode','auto')
[~,b]=fileparts(ni.fname);
[~,b]=fileparts(b);
% print('-painters','-r300','-dpng',[data_path subj_name '/figures/' b '_BBcurves_view' int2str(sl_plot(2)) '_slice' int2str(sl_plot(1))])
print('-painters','-r300','-dpng',[data_path subj_name '/figures/' b '_BBcurves_view' int2str(sl_plot(2)) '_slice' int2str(sl_plot(1)) '_color'])

