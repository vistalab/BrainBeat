clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%

%% Base data directory

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';
% chdir(dDir)

%% Now plot the average responses across subjects

all_subs = [1 2 3 4 5 6];
all_scans = [3 3 3 1 1 1];

[roiCodes,roiNames] = dkt_labels();

t_hr = linspace(-.5,1.5,128);
avResp_hr = NaN(length(t_hr),length(roiNames),length(all_subs));
ppgCurve_hr = NaN(length(t_hr),length(all_subs));
    
for ss = 1:length(all_subs)
    % Functionals
    % Select a subject and scan nummer
    s_nr = all_subs(ss); %[1 2 3 4 5 6]
    scan_nr = all_scans(ss); % [3 3 3 1 1 1]

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

    % get physio stuff we need:
    physio      = physioCreate('nifti',ni);
    ppgCurve    = physioGet(physio,'ppg ppgcurve');
    ppgCurveT   = physioGet(physio,'ppg ppgtcurve');
    srate       = 1/bbGet(ni,'tr');

    ppgResp = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse.nii.gz']));
    ppgT = load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT.mat']));

    niSegm = niftiRead(fullfile(dDir,subj,scan,[scanName '_r_DKTatlas_aseg.nii.gz']));
    % matrix to vector
    segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

    % put the data in a matrix Voxel X Time
    respMat = reshape(ppgResp.data,[size(ppgResp.data,1) * size(ppgResp.data,2) * size(ppgResp.data,3)],size(ppgResp.data,4));

    avResp = zeros(size(ppgResp.data,4),length(roiNames));
    for kk = 1:length(roiNames)
        avResp(:,kk) = mean(respMat(ismember(segmVect,roiCodes{kk}),:),1);
    end

    % resample avResp to ppg cycle
    ppg_cycle   = 1./physioGet(physio,'PPGrate');
    avResp_sel = avResp(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    t_sel = ppgT.t(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle);
    ppgCurve_sel = ppgCurve(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);
    ppgt_sel = ppgCurveT(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);

    t_temp = linspace(min(t_sel),max(t_sel),128);

    avResp_hr(:,:,ss) = interp1(t_sel,avResp_sel,t_temp);
    ppgCurve_hr(:,ss) = interp1(ppgt_sel,ppgCurve_sel,t_temp);

end

%% get distance matric for MDS, try each subject
% X = nanmean(avResp_hr(t_hr<0.5,:,:),3)';

s_nr = 2;
X = avResp_hr(t_hr<0.5,:,s_nr)';

missing_areas = isnan(X(:,1));

% clear missing areas
X(missing_areas,:) = [];
roiUsed = roiNames;
roiUsed(missing_areas) = [];

t_D = t_hr(t_hr<0.5);

D_length = sqrt(sum(X.^2,2));

D_dot = (X*X');

D = D_dot./(D_length*D_length');

D(D>1) = 1;

D_mat = rad2deg(acos(D));

[Y,eigvals] = cmdscale(round(D_mat));

%
figure
subplot(1,2,1)
plot(eigvals)
subplot(1,2,2)
plot(Y(:,1),Y(:,2),'.')
text(Y(:,1),Y(:,2),roiUsed)
% % xlabel('Miles')
% % ylabel('Miles')


%% get distance matric for MDS, try average across subjects
X = nanmean(avResp_hr(t_hr<0.5,:,:),3)';

missing_areas = isnan(X(:,1));

% clear missing areas
X(missing_areas,:) = [];
roiUsed = roiNames;
roiUsed(missing_areas) = [];

t_D = t_hr(t_hr<0.5);

D_length = sqrt(sum(X.^2,2));

D_dot = (X*X');

D = D_dot./(D_length*D_length');

D(D>1) = 1;

D_mat = rad2deg(acos(D));

% [Y,eigvals] = cmdscale(round(D_mat));
[Y] = mdscale(round(D_mat),2); % nonclassical

%
figure
plot(Y(:,1),Y(:,2),'.')
text(Y(:,1),Y(:,2),roiUsed)



%%
areas_plot = [1:4];
figure('Position',[0 0 120 700])
for kk = 1:length(areas_plot)
    thisArea = areas_plot(kk);
    subplot(length(areas_plot)+1,1,kk),hold on
    plot(t_hr,squeeze(avResp_hr(:,thisArea,:))*100,'Color',[0 .6 .8])
    plot(t_hr,mean(avResp_hr(:,thisArea,:),3)*100,'Color',[0 .6 .8],'LineWidth',2)
    plot([0 0],[min(mean(avResp_hr(:,thisArea,:),3)*100)*2 max(mean(avResp_hr(:,thisArea,:),3)*100)*2],'k')
    xlim([min(t_hr) max(t_hr)])
    title(roiNames{thisArea})
end

subplot(length(areas_plot)+1,1,length(areas_plot)+1),hold on
plot(t_hr,ppgCurve_hr./nanstd(ppgCurve_hr,[],1),'Color',[0 .6 .8])
plot(t_hr,mean(ppgCurve_hr./nanstd(ppgCurve_hr,[],1),2),'Color',[0 .6 .8],'LineWidth',2)
plot([0 0],[-2 4],'k')
xlim([min(t_hr) max(t_hr)])
title('PPG')

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir '/figures/segmentation/allsubjects_FA48_tracesDKT1'])
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/allsubjects_FA48_tracesDKT1'])

