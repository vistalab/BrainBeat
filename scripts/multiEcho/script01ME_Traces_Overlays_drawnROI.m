clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';
% chdir(dDir)

%% The T2* data are here.  

% Select a subject and scan nummer
s_nr = 5;
scan_nr = [4 5]; % provide a pair of scans with echo 1 and 2

subs = bb_subs(s_nr);
subj = subs.subj;
% Identify the echo 1 (8ms) and echo 2 (28ms)
scan1 = subs.scan{scan_nr(1)};
scanName1 = subs.scanName{scan_nr(1)};
scan2 = subs.scan{scan_nr(2)};
scanName2 = subs.scanName{scan_nr(2)};

% Load first echo
fmri1 = fullfile(dDir,subj,scan1,[scanName1 '.nii.gz']);
if ~exist(fmri1,'file')
    clear ni1
    error('filename %s does not exist',fmri1)
end
ni1 = niftiRead(fmri1);

% Load second echo
fmri2 = fullfile(dDir,subj,scan2,[scanName2 '.nii.gz']);
if ~exist(fmri2,'file')
    clear ni2
    error('filename %s does not exist',fmri2)
end
ni2 = niftiRead(fmri2);

% Load coregistration matrix for the functionals:
load(fullfile(dDir,subj,scan1,[scanName1 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

% Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,subs.veno,[subs.venoName '.nii']));
% % load coregistration matrix (for the venogram):
% xf_veno=load(fullfile(dDir,subj,subs.veno,[subs.venoName 'AcpcXform.mat']));

%% Anatomicals

anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii.gz']);
niAnatomy = niftiRead(anat);

%% quick overlay between functionals and anatomy

sliceThisDim = 3;
imDims = [-90 -120 -100; 90 130 110];

if s_nr == 5
    curPos = [0,4,38];
end
niFunc = ni1;

% niFunc.data = ni.data(:,:,:,1); % overlay the first functional - more structure visible
% bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
% title('First functional on anatomy')
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/checkCoreg/' subj '_' scan '_Func1onAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

niFunc.data = mean(ni1.data(:,:,:,5:end),4); % overlay the mean functional
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
title('Mean functional on anatomy')
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/checkCoreg/' subj '_' scan '_MeanFuncOnAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])


%% Plot physiology data

figure('Position',[0 0 600 600])

% get data
physio      = physioCreate('nifti',ni1);
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
ppgCurve = physioGet(physio,'ppg ppgcurve');
ppgCurveT = physioGet(physio,'ppg ppgtcurve');
ppgSrate = physioGet(physio,'PPGsrate'); %  = physio.ppg.srate
plot(ppgCurveT,ppgCurve,'k','LineWidth',2);
xlim([ppgCurveT(1) ppgCurveT(end)])
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
respCurve = physioGet(physio,'resp respcurve');
respCurveT = physioGet(physio,'resp resptcurve');
respSrate = physioGet(physio,'RESPsrate');
plot(respCurveT,respCurve,'k','LineWidth',2);
xlabel('time (s)')
xlim([respCurveT(1) respCurveT(end)])
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
% print('-painters','-r300','-depsc',[dDir './figures/physio/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_physioTrace'])
% print('-painters','-r300','-dpng',[dDir './figures/physio/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_physioTrace'])


%% Overlay timeseries

in_data = 'PPG';

sliceThisDim = 1;
imDims = [-90 -120 -100; 90 130 110];
if s_nr == 5 %subject number
    curPos = [6,18,38];
end

% load time series 1 and associated time
ppgTS1 = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_' in_data 'trigResponse.nii.gz'])); % ppg triggered time series
% load time t
load(fullfile(dDir,subj,scan1,[scanName1 '_' in_data 'trigResponseT.mat']),'t');

% load time series 2 and associated time
ppgTS2 = niftiRead(fullfile(dDir,subj,scan2,[scanName2 '_' in_data 'trigResponse.nii.gz'])); % ppg triggered time series
% load time t
load(fullfile(dDir,subj,scan2,[scanName2 '_' in_data 'trigResponseT.mat']),'t');

% load correlation between even and odd scans for colors of timeseries
% ppgRname = fullfile(dDir,subj,scan,[scanName '_corr' in_data '.nii.gz']);
% ppgR = niftiRead(ppgRname); % correlation with PPG
% ppgR.data = ppgR.data.^2;
ppgR1 = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_cod' in_data '.nii.gz'])); % correlation with PPG
ppgR2 = niftiRead(fullfile(dDir,subj,scan2,[scanName2 '_cod' in_data '.nii.gz'])); % correlation with PPG

% %%%% Overlay 1: functionals and anatomy
% niFunc = ni;
% niFunc.data = mean(ni.data(:,:,:,5:end),4); % overlay the mean functional
% bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)

%%%% Overlay 2: timeseries and anatomy
ppgTSplot = ppgTS1;
if isequal(in_data,'PPG')
    ppgTSplot.data(:,:,:,t<-.1 | t>1) = []; % plot these times from curve
elseif isequal(in_data,'RESP')
    ppgTSplot.data(:,:,:,t<-.2 | t>4) = []; % plot these times from curve
end

% Scale factor for size and color:
niColor = ppgR1; % use R2 map to scale later: r2_scale=sqrt(imgSlice2.^2);

% Scale time series amplitude by R, plots are generated with respect to the maximum.
maxTS = max(abs(ppgTSplot.data),[],4); % get the max of each curve
ppgTSplot.data = bsxfun(@rdivide,ppgTSplot.data,maxTS); % devide by the max of each curve (sets all curves to 1 max)
ppgTSplot.data = bsxfun(@times,ppgTSplot.data,niColor.data); % multiply by r^2 to set less reliable curves to zero

bbOverlayTimeseriesAnat(ppgTSplot,niColor,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
% the plotted X and Y do not correspond to actual xyz/mm coordinates, but
% are in a different frame
title('Colors and size scaled by the COD(R)')

clear niColor ppgTSplot
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/voxelTimeSeries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_TraceOnAnat_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])

% bbOverlayTimeseriesVeno(ppgTSplot,niColor,niVeno,acpcXform,xf_veno.acpcXform,sliceThisDim,imDims,curPos)

%% Get ROI timeseries

% Now say we have a voxel let's pull out the raw timeseries and show
% the processing we did:

in_data = 'PPG';

roi_list = {'CFlowvoids','AnteriorSSS','SSS','LeftTransverse','RightTransverse'};

for roi_ind = 1:5%
    niROIname  = roi_list{roi_ind};
    niROI      = niftiRead(fullfile(dDir,subj,subs.anat,['r' subs.anatName niROIname '.nii']));
    
    % get ROI indices:
    [xx,yy,zz] = ind2sub(size(niROI.data),find(niROI.data>0.5));

    % now ROI indices to ACPC (mm):
    xyz_acpc = mrAnatXformCoords(niROI.qto_xyz, [xx,yy,zz]);
    clear xx yy zz % housekeeping

    % now ACPC coordinates to functional indices:
    ijk_func = mrAnatXformCoords(inv(acpcXform), xyz_acpc);
    ijk_func = round(ijk_func); % round to get indices
    ijk_func = unique(ijk_func,'rows'); % only take unique voxels

%     %%%% check for coordinates in functional space
%     z_slice = [5 10 15 18 20 23 25 26 27 30 33 36];
%     figure
%     for kk = 1:length(z_slice)       
%         subplot(3,4,kk)
%         imagesc(ni1.data(:,:,z_slice(kk),1)),hold on
%         xyz_plot = ijk_func(ijk_func(:,3)==z_slice(kk),:);
%         plot(xyz_plot(:,2),xyz_plot(:,1),'r.')
%         axis image 
%     end
    % get ROI curves
    % imgVol = out(1).weights; % PC1 weight
    % imgVol = out(2).weights; % PC2 weight
%     imgVol = svdResults.model;
    imgVol1 = ppgTS1.data;
    imgVol2 = ppgTS2.data;
    % Mask by reliability R
%     imgVol1(ppgR1.data<.5) = NaN;
%     imgVol2(ppgR1.data<.5) = NaN;
    
    roiTrace1 = NaN(size(ijk_func,1),size(imgVol1,4));
    roiTrace2 = NaN(size(ijk_func,1),size(imgVol2,4));
    svdVals = NaN(size(ijk_func,1),2);

    for kk = 1:size(ijk_func,1)
        roiTrace1(kk,:) = imgVol1(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3),:);
        roiTrace2(kk,:) = imgVol2(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3),:);
        % note that x and y are only switched around for plotting position, 
        % not for getting the actual image values
        % svdVals(kk,:) = [out(1).weights(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3)) out(2).weights(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3))];
    end
    % Remove voxels with low R
    roiTrace1(isnan(roiTrace1(:,1)),:) = [];
    roiTrace2(isnan(roiTrace2(:,1)),:) = [];

    figure('Position',[0 0 800 500])
    %%%% Imagesc:
    subplot(2,2,1)
    imagesc(t,[1:size(roiTrace1,1)],roiTrace1,[-1 1])
    subplot(2,2,3)
    imagesc(t,[1:size(roiTrace2,1)],roiTrace2,[-1 1])
    
    %%%% Plot all traces:
%     hold on
%     plot(t,roiTrace')
%     plot([0 0],[min(roiTrace(:)) max(roiTrace(:))],'k')

    %%%% Mesh or surf
%     xx = repmat([1:size(roiTrace,1)],length(t),1)';
%     yy = repmat(t,size(roiTrace,1),1);
%     mesh(xx,yy,roiTrace)
%     view(100,30)

%     xlim([t(1) 1.5])
    
    title(['Traces sub ' int2str(s_nr) '  scan ' int2str(scan_nr) ' roi: ' niROIname])
    subplot(2,2,2)
    plot(t,mean(roiTrace1,1),'k')
    subplot(2,2,4)
    plot(t,mean(roiTrace2,1),'k')

%     title('SVD PC1 and PC2 weights')
%     imagesc(svdVals,[-.02 .02])
%     set(gca,'XTick',[1 2],'XTickLabel',{'PC1','PC2'})
end

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/ROI/' subj '_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_ROIsTest'])

% Get coordinates back in acpc:
xyz_acpc_sparse = mrAnatXformCoords(acpcXform, ijk_func);



%%
    

    % load time series for TE1 and associated time
    ppgTS1 = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_' in_data 'trigResponse.nii.gz'])); % ppg triggered time series
    % load time t
    load(fullfile(dDir,subj,scan1,[scanName1 '_' in_data 'trigResponseT.mat']),'t');
    % load std of time series and associated time
    ppgTSstd1 = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_' in_data 'trigResponse_std.nii.gz'])); % ppg triggered time series
    % get the timeseries for the current voxel
    [voxelTs1] = bbGetVoxelTimeseries(ppgTS1,acpcXform,curPos);
    [voxelTsStd1] = bbGetVoxelTimeseries(ppgTSstd1,acpcXform,curPos);

    % load time series for TE2 and associated time
    ppgTS2 = niftiRead(fullfile(dDir,subj,scan2,[scanName2 '_' in_data 'trigResponse.nii.gz'])); % ppg triggered time series
    % load std of time series and associated time
    ppgTSstd2 = niftiRead(fullfile(dDir,subj,scan2,[scanName2 '_' in_data 'trigResponse_std.nii.gz'])); % ppg triggered time series
    % get the timeseries for the current voxel
    [voxelTs2] = bbGetVoxelTimeseries(ppgTS2,acpcXform,curPos);
    [voxelTsStd2] = bbGetVoxelTimeseries(ppgTSstd2,acpcXform,curPos);

    % plot the timeseries
    figure('Position',[0 0 150 150]),hold on
    plot(t,zeros(size(t)),'Color',[.5 .5 .5])
    plot([0 0],[-2 2],'Color',[.5 .5 .5])

    % % plot TE1 std
    % upErr1 = 100*squeeze(voxelTs1) + 100*voxelTsStd1; % mean + std
    % lowErr1 = 100*squeeze(voxelTs1) - 100*voxelTsStd1;
    % fill([t t(end:-1:1)],[upErr1; lowErr1(end:-1:1)],[.5 .5 1],'EdgeColor',[.5 .5 1],'FaceAlpha',.5)
    % 
    % % plot TE2 std
    % upErr2 = 100*squeeze(voxelTs2) + 100*voxelTsStd2; % mean + std
    % lowErr2 = 100*squeeze(voxelTs2) - 100*voxelTsStd2;
    % fill([t t(end:-1:1)],[upErr2; lowErr2(end:-1:1)],[1 .5 .5],'EdgeColor',[1 .5 .5],'FaceAlpha',.5)

    % plot TE1
    plot(t,100*squeeze(voxelTs1),'b','LineWidth',2)
    % plot TE2
    plot(t,100*squeeze(voxelTs2),'r','LineWidth',2)

    axis tight
    ylabel('% signal modulation') % (signal - mean)./mean
    xlabel('time (s)')
    set(gcf,'PaperPositionMode','auto')
    % print('-painters','-r300','-dpng',[dDir './figures/voxelTimeSeries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_TraceOnAnat_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3)) '_' voxelLabel])
    % print('-painters','-r300','-depsc',[dDir './figures/voxelTimeSeries/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_TraceOnAnat_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3)) '_' voxelLabel])

end

