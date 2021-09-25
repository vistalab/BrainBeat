clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';
% chdir(dDir)

%% Generate basic fMRI change from Freesurfer segmentations 
sub_labels = {'4','5','1'}; 
ses_labels = {'1','1','2'}; 
acq_labels = {'ME','ME','ME'};
run_nrs = {[1],[1],[1]}; % there are 2, 2, 3 runs for these subjects/sessions

dkt_table = readtable('dkt_areas.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dkt_table_surface = readtable('dkt_areas_surface.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
roiNames = dkt_table.label;
roiCodes = dkt_table.label_nr;

t_hr = linspace(-.5,1.5,128);
avResp_hr_lns0 = NaN(length(t_hr),length(roiNames),length(sub_labels));
avRespVeno_hr_lns0 = NaN(length(t_hr),1,length(sub_labels));
avResp_hr_t2s = NaN(length(t_hr),length(roiNames),length(sub_labels));
avRespVeno_hr_t2s = NaN(length(t_hr),1,length(sub_labels));
ppgCurve_hr = NaN(length(t_hr),length(sub_labels));
    
for ss = 1:length(sub_labels)
    % Functionals
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};

    rr = 1;% run_nr
    run_nr = run_nrs{ss}(rr);
    
    % Functionals Echo 1
    fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_echo-1_bold.nii.gz']);
    fmri_name = fullfile(dDir,fmri_BIDSname);
    ni = niftiRead(fmri_name);

    % get physio stuff we need:
    physio      = physioCreate('nifti',ni);
    ppgCurve    = physioGet(physio,'ppg ppgcurve');
    ppgCurveT   = physioGet(physio,'ppg ppgtcurve');
    srate       = 1/bbGet(ni,'tr');
    ppg_cycle   = 1./physioGet(physio,'PPGrate');

    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

    ppgResplnS0 = niftiRead([save_name_base '_PPGtrigResponse_lnS0.nii.gz']); % ppg triggered time series
    ppgRespT2s = niftiRead([save_name_base '_PPGtrigResponse_T2s.nii.gz']); % ppg triggered time series
    ppgT = load([save_name_base '_PPGtrigResponseT.mat'],'t');

    % load segmentation
    niSegm = niftiRead([save_name_base '_r_DKTatlas_aseg.nii.gz']);
    % matrix to vector
    segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

    % also get the Venogram segmentation from niSegm2
    % niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
    niSegm2 = niftiRead([save_name_base '_combineSegm.nii.gz']);
    segm2Vect = reshape(niSegm2.data,[size(niSegm2.data,1) * size(niSegm2.data,2) * size(niSegm2.data,3)],1);
    
    % put the data in a matrix Voxel X Time
    respMat_lnS0 = reshape(ppgResplnS0.data,[size(ppgResplnS0.data,1) * size(ppgResplnS0.data,2) * size(ppgResplnS0.data,3)],size(ppgResplnS0.data,4));
    respMat_t2s = reshape(ppgRespT2s.data,[size(ppgRespT2s.data,1) * size(ppgRespT2s.data,2) * size(ppgRespT2s.data,3)],size(ppgRespT2s.data,4));
    
    % get outliers in T2s and exclude those:
    t2s_outliers = sum(isoutlier(respMat_t2s,'median',1),2)>0;
    
    avResp_lnS0 = zeros(size(ppgResplnS0.data,4),length(roiNames));
    avResp_t2s = zeros(size(ppgRespT2s.data,4),length(roiNames));
    for kk = 1:length(roiNames)
        avResp_lnS0(:,kk) = mean(respMat_lnS0(ismember(segmVect,roiCodes(kk)) & ~t2s_outliers,:),1);
        avResp_t2s(:,kk) = mean(respMat_t2s(ismember(segmVect,roiCodes(kk)) & ~t2s_outliers,:),1);
    end
    
    % get venogram response - remove rows with outliers in median
    avRespVeno_lnS0 = mean(respMat_lnS0(ismember(segm2Vect,5) & ~t2s_outliers,:),1)';
    avRespVeno_t2s = mean(respMat_t2s(ismember(segm2Vect,5) & ~t2s_outliers,:),1)';
    % avRespVeno_t2s = mean(rmoutliers(respMat_t2s(ismember(segm2Vect,5),:),1),1)';

    % resample avResp log s0 to ppg cycle
    avResp_sel_lnS0 = avResp_lnS0(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    avRespVeno_sel_lnS0 = avRespVeno_lnS0(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    % resample avResp t2s to ppg cycle
    avResp_sel_t2s = avResp_t2s(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    avRespVeno_sel_t2s = avRespVeno_t2s(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    % get timing
    t_sel = ppgT.t(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle);
    % resample ppg curve
    ppgCurve_sel = ppgCurve(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);
    ppgt_sel = ppgCurveT(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);

    % resampling to hartbeat time with 128 timepoints t_hr
    t_temp = linspace(min(t_sel),max(t_sel),128);

    avResp_hr_lns0(:,:,ss) = interp1(t_sel,avResp_sel_lnS0,t_temp);
    avRespVeno_hr_lns0(:,:,ss) = interp1(t_sel,avRespVeno_sel_lnS0,t_temp);
    avResp_hr_t2s(:,:,ss) = interp1(t_sel,avResp_sel_t2s,t_temp);
    avRespVeno_hr_t2s(:,:,ss) = interp1(t_sel,avRespVeno_sel_t2s,t_temp);
    ppgCurve_hr(:,ss) = interp1(ppgt_sel,ppgCurve_sel,t_temp);
    
    clear t_sel ppgt_sel t_temp % t_sel and ppgt_sel are subject specific
    clear avResp_sel_lnS0 avResp_sel_t2s avRespVeno_sel_lnS0 avRespVeno_sel_t2s...
        ppgCurve_sel avRespVeno_lnS0 avRespVeno_t2s avResp_lnS0 avResp_t2s...
        ppgCurveT ppgCurve ppg_cycle ppg_onsets % all other subject specific vars
end

%%

%%
%% Figure 4D: plot set of responses

n_subs = size(avResp_hr_lns0,3);

figure('Position',[0 0 300 400])

subplot(4,2,1),hold on
this_area = 1; % caudal anterior cingulate 1 
thisSignal = avResp_hr_lns0(:,this_area,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,2,2),hold on
thisSignal = avResp_hr_t2s(:,this_area,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])
% title(roiNames{this_area});

subplot(4,2,3),hold on
% title('Lateral ventricles')
% average left and lateral ventricles 74 89
thisSignal = (squeeze(avResp_hr_lns0(:,74,:))+squeeze(avResp_hr_lns0(:,89,:)))/2;
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_hr,zscore(thisSignal));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,2,4),hold on
thisSignal = (squeeze(avResp_hr_t2s(:,74,:))+squeeze(avResp_hr_t2s(:,89,:)))/2;
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_hr,zscore(thisSignal));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])


subplot(4,2,5),hold on
% title('Veins')
thisSignal = avRespVeno_hr_lns0(:,1,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,2,6),hold on
% title('Veins')
thisSignal = avRespVeno_hr_t2s(:,1,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

%%
subplot(4,1,4),hold on
% title('PPG')
thisSignal = ppgCurve_hr;
thisSignal(isnan(thisSignal)) = 0;
thisSignal = zscore(squeeze(thisSignal));
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
fill([t_hr t_hr(end:-1:1)],100*[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot([0 0],[-1 3],'Color',[.7 .7 .7])
plot(t_hr,100*thisSignal,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

