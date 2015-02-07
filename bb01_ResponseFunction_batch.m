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

%% Preprocess the data: PPG time series and correlation

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

for s=2:3
    s_info = bb_subs(s);
    subj=s_info.subj;
    for scan_nr=1:length(s_info.scan)
        scan=s_info.scan{scan_nr};
        scanName=s_info.scanName{scan_nr};

        fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
        ni = niftiRead(fmri);

        % compute the PPG triggered response matrix for the whole brain and save as a nifti:
        [response_matrix,t,response_matrix_odd,response_matrix_even] = bbResponse2physio(ni);

        % safe time T:
        save(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT']),'t')

        % save average of all heartbeats:
        ni1=ni;
        ni1.data=response_matrix;
        ni1.fname=[scanName '_PPGtrigResponse'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
        
        % save average of all odd heartbeats:
        ni1=ni;
        ni1.data=response_matrix_odd;
        ni1.fname=[scanName '_PPGtrigResponse_odd'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
        
        % save average of all even heartbeats:
        ni1=ni;
        ni1.data=response_matrix_even;
        ni1.fname=[scanName '_PPGtrigResponse_even'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1

    end
    for scan_nr=1:length(s_info.scan)
        scan=s_info.scan{scan_nr};
        scanName=s_info.scanName{scan_nr};

        ni_odd = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_odd']));
        ni_even = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_even']));
        
        % compute the correlation/ reliability with the PPG for the whole brain and save as a nifti:
        [out_r_map]=bbCorrelate2physio(ni_odd,ni_even);

        % save as nifti
        ni1=ni;
        ni1.data=out_r_map;
        ni1.fname=[scanName '_corrPPG.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
    end
end
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

