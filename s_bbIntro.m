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
scan ='6_1_mux8fov4_r1_25s_4mmFA25';  % A particular data set
fmri = fullfile(dDir,subj,scan,'8202_6_1.nii.gz');
ni = niftiRead(fmri);


%% We should get the anatomicals up too at some point

anat      ='9_1_T1w_1mm_sag';   % Anatomical data
anat      = fullfile(dDir,subj,anat,'8202_9_1.nii.gz');
niAnatomy = niftiRead(anat);

% mrViewer should be able to take a nifti structure as input data
% The orientation of the 3-axis images is in Radiological format.  Perhaps
% this should be the default?
% The zoom is incomprehensible or broken.
mrViewer(anat,'nifti')

%% Let's have a look at some of the key parameters

BB_get(ni,'timing')



%% Let's pick a physiology file and do something

% load the physiology data:
physio = BB_get(ni,'physio');

% let's get the ppg peaks in seconds, this also produces a plot of the ppg
% data locked to the peak
ppg_onsets=BB_get(ni,'ppg_peaks');

%% Let's look at the correlation with PPG in a slice

% pick a slice
slice_plot = 20;

% calculate the correlation
out_r_map = BB_get(ni,'ppg_correlation',slice_plot);

% make a figure
figure('Position',[0 0 800 300])
subplot(1,2,1), imagesc(ni.data(:,:,sl_plot,1))
axis image, title('volume 1, for reference')
subplot(1,2,2), imagesc(out_r_map,[0 1])
colorbar, axis image, title('correlation with heartbeat (r)')

%% Let's make a movie of the signal after PPG peak

sl_plot=20;

% calculate the response matrix
val = BB_get(ni,'ppg_response_function',sl_plot);
response_matrix = val.response_matrix;
t = val.t;
clear val

% get ppg response curve with same timing as response_matrix for plotting
ppg_ind = 1; % PPG 
phys_srate = physio(ppg_ind).CNIsrate;
signal = physio(ppg_ind).data;
epoch_pre = -t(1)*phys_srate;
epoch_post = t(end)*phys_srate;
ppg_t = (-round(epoch_pre)+1:epoch_post)/phys_srate;

% get rid of early and late ppg onsets that do have an entire curve:
timing = BB_get(ni,'timing');% we need the timing
ppg_onsets=ppg_onsets((ppg_onsets+t(1))>0); % get rid of early ones
ppg_onsets=ppg_onsets((ppg_onsets+t(end))<max(timing(:))); % get rid of late ones

phys_epochs1 = zeros(length(ppg_onsets),length(ppg_t)); % create epochs of each ppg response
for k=1:length(ppg_onsets)
    sample_ppg=round(ppg_onsets(k)*phys_srate);
    phys_epochs1(k,:)=signal(sample_ppg-round(epoch_pre)+1:sample_ppg+(epoch_post));
end
ppg_resp=mean(phys_epochs1,1);% get average ppg response
clear phys_epochs sample_ppg ppg_onsets signal % housekeeping

% initiate variables for the movie
nFrames=length(t);
mov(1:nFrames)=struct('cdata',[],'colormap',[]);

% make the figure for the movie
figure('Position',[0 0 400 400])
subplot(2,1,2),hold on
plot(ppg_t,ppg_resp-mean(ppg_resp)); 
xlabel('time (s)'),ylabel('PPG response')

% subtract the mean, to look at changes over the variance
response_matrix_meansub=response_matrix-repmat(mean(response_matrix,3),[1 1 size(response_matrix,3)]);

for k=1:size(response_matrix,3)
    subplot(2,2,1)
    imagesc(response_matrix(:,:,k),[0 max(response_matrix(:))])
    axis image, title(['data t=' num2str(t(k))])
    
    subplot(2,2,2)
    imagesc(response_matrix_meansub(:,:,k),[-2 2])
    axis image, title(['data-mean t=' num2str(t(k))])

    subplot(2,1,2),hold on
    xlim([t(1) t(end)]), plot(t(k),0,'r.')
    
    pause(1)
    mov(k)=getframe(gcf);
end
close gcf

% [~,b]=fileparts(ni.fname);
% [~,b]=fileparts(b);
movie2avi(mov,fullfile(dDir,subj,scan,[b '_PPGmov_slice' int2str(sl_plot) '.avi']),'compression','None','fps',1)


%% End