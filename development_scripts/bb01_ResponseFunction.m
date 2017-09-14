clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%
%

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
dDir = '/biac4/wandell/data/BrainBeat/data';

% chdir(dDir)

%% The T2* data are here.  

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

subj ='20141017_1242';    % Date _ Time out of NIMS
scan ='5_1_mux8fov4_r1_25s_4mm';  % A particular data set
scanName='8202_5_1';
fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
ni = niftiRead(fmri);

% subj_name='20141017_1242';
% scan_name='5_1_mux8fov4_r1_25s_4mm';
% f_name='8202_5_1';
% scan_name='6_1_mux8fov4_r1_25s_4mmFA25';
% f_name='8202_6_1';
% scan_name='7_1_mux8fov4_r1_25s_4mmFA48';
% f_name='8202_7_1';

% subj_name='20140911_1411';
% scan_name='3_1_mux8fov4_r1_35s_3mm';
% f_name='7944_3_1';
% scan_name='6_1_mux8fov4_r1_25s_4mm';
% f_name='7944_6_1';

%% We should get the anatomicals up too at some point

anat      ='9_1_T1w_1mm_sag';   % Anatomical data
anat      = fullfile(dDir,subj,anat,'8202_9_1.nii.gz');
niAnatomy = niftiRead(anat);

%% plot physiology data

physio     = physioCreate('nifti',ni,'figure',1);
subplot(2,1,1),xlim([0 20])
subplot(2,1,2),xlim([0 20])
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/2014_12_presentation1/' subj '_' scan '_physioraw'])


ppgData = physioGet(physio,'ppg data');
respData = physioGet(physio,'resp data');

% get PPG peaks
ppgPeaks = physioGet(physio,'ppg peaks');

% plot PPG peaks
t = physioGet(physio,'ppg sample times');
peakVal = max(ppgData(:))*ones(size(ppgPeaks));
mrvNewGraphWin;
srate = physioGet(physio,'ppg srate');
peakSamples = round(ppgPeaks*srate);
plot(ppgPeaks,ppgData(peakSamples),'ro',t,ppgData,'k-');
xlabel('secs'); grid on
% xlim([0 20])
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/2014_12_presentation1/' subj '_' scan '_ppgPeaks'])

% get RESP peaks
respPeaks = physioGet(physio,'resp peaks');

% plot RESP peaks
t = physioGet(physio,'resp sample times');
peakVal = max(respData(:))*ones(size(respPeaks));
mrvNewGraphWin;
srate = physioGet(physio,'resp srate');
peakSamples = round(respPeaks*srate);
plot(respPeaks,respData(peakSamples),'ro',t,respData,'k-');
xlabel('secs'); grid on
% xlim([0 40])
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/2014_12_presentation1/' subj '_' scan '_respPeaks'])

% get PPG rate
ppgRate = physioGet(physio,'ppg rate');

% get RESP rate
respRate = physioGet(physio,'resp rate');

% get one PPG curve
ppgCurve = physioGet(physio,'ppg curve');

%% get responses after PPG peak

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
%% compute the correlation/ reliability with the PPG for the whole brain and save as a nifti:
%%

[out_r_map]=bbCorrelate2physio(ni);

ni1=ni;
ni1.data=out_r_map;
ni1.fname=[scanName '_corrPPG.nii.gz'];
niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
clear ni1

%%
%% compute the PPG triggered response matrix for the whole brain and save as a nifti:
%%

% now run for the whole brain:
[response_matrix,t] = bbResponse2physio(ni);

% safe the response matrix as nifti and save  time t
ni1=ni;
ni1.data=response_matrix;
ni1.fname=[scanName '_PPGtrigResponse'];
niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
clear ni1
save(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT']),'t')


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
