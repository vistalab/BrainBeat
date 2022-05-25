
% code to make Figure 2

clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%%
%% Overlay anat and reliability (COD)
%%

sub_labels = {'1'}; 
ses_labels = {'1'}; 
acq_labels = {'4mmFA48'};
run_nrs = {[1]};

ss = 1;%:length(sub_labels) % subjects/ses/acq
rr = 1;% run_nr
sub_label = sub_labels{ss};
ses_label = ses_labels{ss};
acq_label = acq_labels{ss};

run_nr = run_nrs{ss}(rr);
        
% Get the anatomicals:
t1w_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'anat',...
            ['sub-' sub_label '_ses-' ses_label '_T1w.nii.gz']);
niAnatomy = niftiRead(fullfile(dDir,t1w_BIDSname));

save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);

% load coregistration matrix (for the functionals):
load([save_name_base '_AcpcXform_new.mat']);
acpcXform = acpcXform_new; clear acpcXform_new

%%%% Overlay 1: cod and anatomy
ppgRname = [save_name_base '_codPPG.nii.gz'];
ppgR = niftiRead(ppgRname); % correlation with PPG

%% Subject 1 Sagittal slices T1 + COD
curPos = [-10,1,-20]; 
sliceThisDim = 1; 
imDims = [-90 -120 -120; 90 130 90];
overlayPlot = ppgR;
cod_th = 0.5;

bb_roi = bb_subs_rois(ss); % script that lists voxel indices for subject 1

for kk = -1%-60:10:60
    curPos(1) = kk;
    bbOverlayFuncAnat(overlayPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,1,1,cod_th);
    title(['slice ' int2str(kk) ', R^2>' num2str(cod_th,3)])

    % add positive for kk = 1
    if kk == -1
        % Loop through ROI indices/labels
        for roi_ind = 1:8
            curPos = bb_roi(roi_ind).curPos;
            voxelLabel = bb_roi(roi_ind).voxelLabel;
            plot(curPos(2),curPos(3),'w.')
        end
    end
    set(gcf,'PaperPositionMode','auto')    
    print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
        ['Figure2A_sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_orient' int2str(sliceThisDim) '_slice' int2str(kk)]))
    print('-r300','-depsc2',fullfile(dDir,'derivatives','figures',...
        ['Figure2A_sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_orient' int2str(sliceThisDim) '_slice' int2str(kk)]))
end

%% Plot timeseries

% load COD
ppgRname = [save_name_base '_codPPG.nii.gz'];

% ppg locked time series
ppgTS = niftiRead([save_name_base '_PPGtrigResponse.nii.gz']); 
% associated time t
load([save_name_base '_PPGtrigResponseT.mat'],'t');
% load std of time series
ppgTSsterr = niftiRead([save_name_base '_PPGtrigResponse_sterr.nii.gz']); % ppg triggered time series

bb_roi = bb_subs_rois(ss); % script that lists voxel indices

for roi_ind = 1:8
    % Call ROI indices/labels
    curPos = bb_roi(roi_ind).curPos;
    voxelLabel = bb_roi(roi_ind).voxelLabel;
    % get the timeseries
    [voxelTs] = bbGetVoxelTimeseries(ppgTS,acpcXform,curPos);
    [voxelTsSterr] = bbGetVoxelTimeseries(ppgTSsterr,acpcXform,curPos);
    voxelTsSterr = 2*voxelTsSterr; % 2*standard error

    % plot the timeseries
    figure('Position',[0 0 130 120]),hold on
    plot(t,zeros(size(t)),'Color',[.5 .5 .5])
    plot([0 0],[-2 2],'Color',[.5 .5 .5])
    upErr = 100*squeeze(voxelTs) + 100*voxelTsSterr; % mean + 2*standard error
    lowErr = 100*squeeze(voxelTs) - 100*voxelTsSterr;
    fill([t t(end:-1:1)],[upErr; lowErr(end:-1:1)],[.5 .5 .5],'EdgeColor',[.5 .5 .5])
    plot(t,100*squeeze(voxelTs),'k','LineWidth',2)
    axis tight
    ylabel('% signal modulation') % (signal - mean)./mean
    xlabel('time (s)')
    title(bb_roi(roi_ind).voxelLabel)
    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
        ['Figure2_sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3)) '_' voxelLabel]))
    print('-r300','-depsc2',fullfile(dDir,'derivatives','figures',...
        ['Figure2_sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_PosMM' int2str(curPos(1)) '_' int2str(curPos(2)) '_' int2str(curPos(3)) '_' voxelLabel]))
end

%%
%%
%% %% Overlay 2: timeseries and anatomy (not used)
%%
%%
% 
% ppgTSplot = ppgTS;
% ppgTSplot.data(:,:,:,t<-.1 | t>1)=[]; % plot these times from curve
% 
% % Scale time series amplitude by R, plots are generated with respect to the maximum.
% niColor = ppgR; % use R2 map to scale later: r2_scale=sqrt(imgSlice2.^2);
% maxTS = max(abs(ppgTSplot.data),[],4); % get the max of each curve
% ppgTSplot.data = bsxfun(@rdivide,ppgTSplot.data,maxTS); % devide by the max of each curve (sets all curves to 1 max)
% ppgTSplot.data = bsxfun(@times,ppgTSplot.data,niColor.data); % multiply by r^2 to set less reliable curves to zero
% 
% bbOverlayTimeseriesAnat(ppgTSplot,niColor,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)
% 
% clear niColor ppgTSplot
% 
