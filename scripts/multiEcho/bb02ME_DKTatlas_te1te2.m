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
run_nrs = {[2],[2],[2]}; % there are 2, 2, 3 runs for these subjects/sessions

dkt_table = readtable('dkt_areas.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dkt_table_surface = readtable('dkt_areas_surface.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
roiNames = dkt_table.label;
roiCodes = dkt_table.label_nr;

t_hr = linspace(-.5,1.5,128);
avResp_hr_TE1 = NaN(length(t_hr),length(roiNames),length(sub_labels));
avRespVeno_hr_TE1 = NaN(length(t_hr),1,length(sub_labels));
avResp_hr_TE2 = NaN(length(t_hr),length(roiNames),length(sub_labels));
avRespVeno_hr_TE2 = NaN(length(t_hr),1,length(sub_labels));
avResp_hr_S0 = NaN(length(t_hr),length(roiNames),length(sub_labels));
avRespVeno_hr_S0 = NaN(length(t_hr),1,length(sub_labels));
avResp_hr_T2s = NaN(length(t_hr),length(roiNames),length(sub_labels));
avRespVeno_hr_T2s = NaN(length(t_hr),1,length(sub_labels));
ppgCurve_hr = NaN(length(t_hr),length(sub_labels));
    
for ss = 1:3%:length(sub_labels)
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
    ni1 = niftiRead(fmri_name);
    % Functionals Echo 2
    fmri_BIDSname2 = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_echo-2_bold.nii.gz']);
    fmri_name2 = fullfile(dDir,fmri_BIDSname2);
    ni2 = niftiRead(fmri_name2);

    % get physio stuff we need:
    physio      = physioCreate('nifti',ni1);
    ppgCurve    = physioGet(physio,'ppg ppgcurve');
    ppgCurveT   = physioGet(physio,'ppg ppgtcurve');
    srate       = 1/bbGet(ni1,'tr');
    ppg_cycle   = 1./physioGet(physio,'PPGrate');

    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)]);

    ppgRespTE1 = niftiRead([save_name_base '_PPGtrigResponse_TE1.nii.gz']); % ppg triggered time series
    ppgRespTE2 = niftiRead([save_name_base '_PPGtrigResponse_TE2.nii.gz']); % ppg triggered time series
    ppgT = load([save_name_base '_PPGtrigResponseTE12.mat'],'t');

    % load segmentation
    niSegm = niftiRead([save_name_base '_r_DKTatlas_aseg.nii.gz']);
    % matrix to vector
    segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

    % also get the Venogram segmentation from niSegm2
    % niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
    niSegm2 = niftiRead([save_name_base '_combineSegm.nii.gz']);
    segm2Vect = reshape(niSegm2.data,[size(niSegm2.data,1) * size(niSegm2.data,2) * size(niSegm2.data,3)],1);
    
    % put the data in a matrix Voxel X Time
    respMat_TE1 = reshape(ppgRespTE1.data,[size(ppgRespTE1.data,1) * size(ppgRespTE1.data,2) * size(ppgRespTE1.data,3)],size(ppgRespTE1.data,4));
    respMat_TE2 = reshape(ppgRespTE2.data,[size(ppgRespTE2.data,1) * size(ppgRespTE2.data,2) * size(ppgRespTE2.data,3)],size(ppgRespTE2.data,4));
    
    % now we can calculate S0 and T2s here
    % get echo times in units of seconds
    te1 = 0.001*bbGet(ni1,'te');
    te2 = 0.001*bbGet(ni2,'te');
    % calculate T2*
    rr = (respMat_TE2 - respMat_TE1)/(te2 - te1);
    t2s = -1./rr;
    % calculate log S0, exponentiate later?
    lns0 = respMat_TE1 - (respMat_TE2 - respMat_TE1)*(te1/(te2 - te1));
    
    % get outliers in TE2 and exclude those:
    te1_outliers = sum(isoutlier(respMat_TE1,'median',1),2)>0;
    
    % initiate averages per ROI 
    avResp_TE1 = NaN(size(ppgRespTE1.data,4),length(roiNames));
    avResp_TE2 = NaN(size(ppgRespTE2.data,4),length(roiNames));
    avResp_S0 = NaN(size(ppgRespTE1.data,4),length(roiNames));
    avResp_T2s = NaN(size(ppgRespTE2.data,4),length(roiNames));
    for kk = 1:length(roiNames)
%         avResp_TE1(:,kk) = mean(respMat_TE1(ismember(segmVect,roiCodes(kk)) & ~te1_outliers,:),1);
%         avResp_TE2(:,kk) = mean(respMat_TE2(ismember(segmVect,roiCodes(kk)) & ~te1_outliers,:),1);
        avResp_TE1(:,kk) = mean(respMat_TE1(ismember(segmVect,roiCodes(kk)),:),1);
        avResp_TE2(:,kk) = mean(respMat_TE2(ismember(segmVect,roiCodes(kk)),:),1);
%         avResp_S0(:,kk) = mean(lns0(ismember(segmVect,roiCodes(kk)) & ~te1_outliers,:),1);
%         avResp_T2s(:,kk) = mean(t2s(ismember(segmVect,roiCodes(kk)) & ~te1_outliers,:),1);
        avResp_S0(:,kk) = mean(lns0(ismember(segmVect,roiCodes(kk)),:),1);
        avResp_T2s(:,kk) = mean(t2s(ismember(segmVect,roiCodes(kk)),:),1);
    end
    
    % get venogram response - remove rows with outliers in median
%     avRespVeno_TE1 = mean(respMat_TE1(ismember(segm2Vect,5) & ~te1_outliers,:),1)';
%     avRespVeno_TE2 = mean(respMat_TE2(ismember(segm2Vect,5) & ~te1_outliers,:),1)';
%     avRespVeno_S0 = mean(lns0(ismember(segm2Vect,5) & ~te1_outliers,:),1)';
%     avRespVeno_T2s = mean(t2s(ismember(segm2Vect,5) & ~te1_outliers,:),1)';
    avRespVeno_TE1 = mean(respMat_TE1(ismember(segm2Vect,5),:),1)';
    avRespVeno_TE2 = mean(respMat_TE2(ismember(segm2Vect,5),:),1)';
    avRespVeno_S0 = mean(lns0(ismember(segm2Vect,5),:),1)';
    avRespVeno_T2s = mean(t2s(ismember(segm2Vect,5),:),1)';
    % avRespVeno_t2s = mean(rmoutliers(respMat_t2s(ismember(segm2Vect,5),:),1),1)';

    % resample avResps to ppg cycle
    avResp_sel_TE1 = avResp_TE1(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    avRespVeno_sel_TE1 = avRespVeno_TE1(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    avResp_sel_TE2 = avResp_TE2(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    avRespVeno_sel_TE2 = avRespVeno_TE2(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    avResp_sel_S0 = avResp_S0(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    avRespVeno_sel_S0 = avRespVeno_S0(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    avResp_sel_T2s = avResp_T2s(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    avRespVeno_sel_T2s = avRespVeno_T2s(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);

    % get timing
    t_sel = ppgT.t(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle);
    % resample ppg curve
    ppgCurve_sel = ppgCurve(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);
    ppgt_sel = ppgCurveT(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);

    % resampling to hartbeat time with 128 timepoints t_hr
    t_temp = linspace(min(t_sel),max(t_sel),128);

    avResp_hr_TE1(:,:,ss) = interp1(t_sel,avResp_sel_TE1,t_temp);
    avRespVeno_hr_TE1(:,:,ss) = interp1(t_sel,avRespVeno_sel_TE1,t_temp);
    avResp_hr_TE2(:,:,ss) = interp1(t_sel,avResp_sel_TE2,t_temp);
    avRespVeno_hr_TE2(:,:,ss) = interp1(t_sel,avRespVeno_sel_TE2,t_temp);
    avResp_hr_S0(:,:,ss) = interp1(t_sel,avResp_sel_S0,t_temp);
    avRespVeno_hr_S0(:,:,ss) = interp1(t_sel,avRespVeno_sel_S0,t_temp);
    avResp_hr_T2s(:,:,ss) = interp1(t_sel,avResp_sel_T2s,t_temp);
    avRespVeno_hr_T2s(:,:,ss) = interp1(t_sel,avRespVeno_sel_T2s,t_temp);
    ppgCurve_hr(:,ss) = interp1(ppgt_sel,ppgCurve_sel,t_temp);
    
    clear t_sel ppgt_sel t_temp % t_sel and ppgt_sel are subject specific
    clear avResp_sel_TE1 avResp_sel_TE2 avRespVeno_sel_TE1 avRespVeno_sel_TE2...
        avRespVeno_TE1 avRespVeno_TE2 avResp_TE1 avResp_TE2...
        avResp_sel_S0 avResp_sel_T2s avRespVeno_sel_S0 avRespVeno_sel_T2s...
        avRespVeno_S0 avRespVeno_T2s avResp_S0 avResp_T2s...
        ppgCurve_sel ppgCurveT ppgCurve ppg_cycle ppg_onsets % all other subject specific vars
end

%%

%%
%% Figure 4D: plot set of responses TE1 TE2

n_subs = size(avResp_hr_TE1,3);

figure('Position',[0 0 600 400])

subplot(4,2,1),hold on
this_area = 1; % caudal anterior cingulate 1 
thisSignal = avResp_hr_TE1(:,this_area,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])
ylabel(roiNames{this_area}([1 5 end-8:end]));

subplot(4,2,2),hold on
thisSignal = avResp_hr_TE2(:,this_area,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])


subplot(4,2,3),hold on
ylabel('Lateral ventricles')
% average left and lateral ventricles 74 89
thisSignal = (squeeze(avResp_hr_TE1(:,74,:))+squeeze(avResp_hr_TE1(:,89,:)))/2;
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(thisSignal));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(thisSignal));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,2,4),hold on
thisSignal = (squeeze(avResp_hr_TE2(:,74,:))+squeeze(avResp_hr_TE2(:,89,:)))/2;
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(thisSignal));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(thisSignal));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])


subplot(4,2,5),hold on
ylabel('Veins')
thisSignal = avRespVeno_hr_TE1(:,1,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,2,6),hold on
% title('Veins')
thisSignal = avRespVeno_hr_TE2(:,1,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])
for kk = 1:6
    subplot(4,2,kk),hold on
    plot([0 0],[-2 2],'k')
    if mod(kk,2)==1
        title('TE1')
    else
        title('TE2')
    end
end


for ll = 7:8
    subplot(4,2,ll),hold on
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
end

%%

%% Figure 4D: plot set of responses ln S0, T2s

n_subs = size(avResp_hr_TE1,3);

figure('Position',[0 0 600 400])

subplot(4,2,1),hold on
this_area = 1; % caudal anterior cingulate 1 
thisSignal = avResp_hr_S0(:,this_area,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])
ylabel(roiNames{this_area}([1 5 end-8:end]));

subplot(4,2,2),hold on
thisSignal = avResp_hr_T2s(:,this_area,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])


subplot(4,2,3),hold on
ylabel('Lateral ventricles')
% average left and lateral ventricles 74 89
thisSignal = (squeeze(avResp_hr_S0(:,74,:))+squeeze(avResp_hr_S0(:,89,:)))/2;
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(thisSignal));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(thisSignal));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,2,4),hold on
thisSignal = (squeeze(avResp_hr_T2s(:,74,:))+squeeze(avResp_hr_T2s(:,89,:)))/2;
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(thisSignal));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(thisSignal));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])


subplot(4,2,5),hold on
ylabel('Veins')
thisSignal = avRespVeno_hr_S0(:,1,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,2,6),hold on
thisSignal = avRespVeno_hr_T2s(:,1,:);
up_ci = mean(thisSignal,3) + 2*std(thisSignal,[],3)/sqrt(n_subs);
low_ci = mean(thisSignal,3) - 2*std(thisSignal,[],3)/sqrt(n_subs);
% fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% plot(t_hr,(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
plot(t_hr,zscore(squeeze(thisSignal)));%,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

for kk = 1:6
    subplot(4,2,kk),hold on
    plot([0 0],[-2 2],'k')
    if mod(kk,2)==1
        title('S0')
    else
        title('T2s')
    end
end


for ll = 7:8
    subplot(4,2,ll),hold on
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
end