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

dkt_table = readtable(fullfile(dDir,'dkt_areas.tsv'),'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
roiNames = dkt_table.label;
roiCodes = dkt_table.label_nr;

t_hr = linspace(-.5,1.5,128);
avResp_hr = NaN(length(t_hr),length(roiNames),length(all_subs));
avRespVeno_hr = NaN(length(t_hr),1,length(all_subs));
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

    %%%% Scale the time-series matrix by the reliability
    % Get odd/even cod (made with bbCod/Correlate2physio):
    ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
    ppgR = niftiRead(ppgRname); % COD between even and odd heartbeats
    % Set maximum of ppgTS to 1 for each voxel
    ppgResp.data = ppgResp.data ./ repmat(max(abs(ppgResp.data),[],4),[1,1,1,size(ppgResp.data,4)]);
    % Multiply by cod size (always >= 0)
    ppgResp.data = ppgResp.data .* repmat(ppgR.data,[1,1,1,size(ppgResp.data,4)]);
    
    % load segmentation
    niSegm = niftiRead(fullfile(dDir,subj,scan,[scanName '_r_DKTatlas_aseg.nii.gz']));
    % matrix to vector
    segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

    % also get the Venogram segmentation from niSegm2
    % niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
    niSegm2 = niftiRead(fullfile(dDir,subj,scan,[scanName '_combineSegm.nii.gz']));
    segm2Vect = reshape(niSegm2.data,[size(niSegm2.data,1) * size(niSegm2.data,2) * size(niSegm2.data,3)],1);
    
    % put the data in a matrix Voxel X Time
    respMat = reshape(ppgResp.data,[size(ppgResp.data,1) * size(ppgResp.data,2) * size(ppgResp.data,3)],size(ppgResp.data,4));

    avResp = zeros(size(ppgResp.data,4),length(roiNames));
    for kk = 1:length(roiNames)
        avResp(:,kk) = mean(respMat(ismember(segmVect,roiCodes(kk)),:),1);
    end
    
    % get venogram response
    avRespVeno = mean(respMat(ismember(segm2Vect,5),:),1)';

    % resample avResp to ppg cycle
    ppg_cycle   = 1./physioGet(physio,'PPGrate');
    avResp_sel = avResp(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    avRespVeno_sel = avRespVeno(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle,:);
    t_sel = ppgT.t(ppgT.t>=(0-(.5*ppg_cycle)) & ppgT.t<=1.5*ppg_cycle);
    ppgCurve_sel = ppgCurve(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);
    ppgt_sel = ppgCurveT(ppgCurveT>=(0-(.5*ppg_cycle)) & ppgCurveT<=1.5*ppg_cycle);

    t_temp = linspace(min(t_sel),max(t_sel),128);

    avResp_hr(:,:,ss) = interp1(t_sel,avResp_sel,t_temp);
    avRespVeno_hr(:,:,ss) = interp1(t_sel,avRespVeno_sel,t_temp);
    ppgCurve_hr(:,ss) = interp1(ppgt_sel,ppgCurve_sel,t_temp);
    
end

%%
%% plot some responses for Figure 1

n_subs = size(avResp_hr,3);

figure('Position',[0 0 150 280])

subplot(4,1,1),hold on
this_area = 1; % caudal anterior cingulate 1 
% title(roiNames{this_area});
up_ci = mean(avResp_hr(:,this_area,:),3) + 2*std(avResp_hr(:,this_area,:),[],3)/sqrt(n_subs);
low_ci = mean(avResp_hr(:,this_area,:),3) - 2*std(avResp_hr(:,this_area,:),[],3)/sqrt(n_subs);
fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,1,2),hold on
% title('Lateral ventricles')
thisSignal = (squeeze(avResp_hr(:,74,:))+squeeze(avResp_hr(:,89,:)))/2;
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% average left and lateral ventricles 74 89
plot(t_hr,thisSignal,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,1,3),hold on
% title('Veins')
up_ci = mean(avRespVeno_hr(:,1,:),3) + 2*std(avRespVeno_hr(:,1,:),[],3)/sqrt(n_subs);
low_ci = mean(avRespVeno_hr(:,1,:),3) - 2*std(avRespVeno_hr(:,1,:),[],3)/sqrt(n_subs);
fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_hr,squeeze(avRespVeno_hr(:,1,:)),'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,1,4),hold on
% title('PPG')
thisSignal = ppgCurve_hr;
thisSignal(isnan(thisSignal)) = 0;
thisSignal = zscore(squeeze(thisSignal));
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
fill([t_hr t_hr(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot([0 0],[-1 3],'Color',[.7 .7 .7])
plot(t_hr,thisSignal,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir '/figures/segmentation/Fig1D_allSubs_tracessegm01'])
print('-painters','-r300','-depsc',[dDir '/figures/segmentation/Fig1D_allSubs_tracessegm01'])



%% get distance matric for MDS, try individual subjects
% X = nanmean(avResp_hr(t_hr<0.5,:,:),3)';

s_nr = 2;
t_select = find(t_hr<0.5);

X = avResp_hr(t_select,:,s_nr)';
t_D = t_hr(t_select);

% clear missing areas
missing_areas = isnan(X(:,1));
X(missing_areas,:) = [];
roiUsed = roiNames;
roiUsed(missing_areas) = [];


% construct distance matrix
D_length = sqrt(sum(X.^2,2));

D_dot = (X*X');

D = D_dot./(D_length*D_length');

D(D>1) = 1; 

D_mat = rad2deg(acos(D));

% multidimensional scaling
% [Y,eigvals] = cmdscale(round(D_mat));
[Y,eigvals] = mdscale(round(D_mat),2);

%
figure
% subplot(1,2,1)
% plot(eigvals)
% subplot(1,2,2)
plot(Y(:,1),Y(:,2),'.')
text(Y(:,1),Y(:,2),roiUsed)
% % xlabel('Miles')
% % ylabel('Miles')

figure('Position',[0 0 400 1000])
imagesc(t_D,[1:length(roiUsed)],X)
set(gca,'YTick',[1:length(roiUsed)],'YTickLabel',roiUsed,'FontSize',9,'FontName','Arial')


%% get distance matric for MDS, try average across subjects

t_select = find(t_hr<0.5);

X = nanmean(avResp_hr(t_select,:,:),3)';
t_D = t_hr(t_select);


% clear missing areas
missing_areas = isnan(X(:,1));
X(missing_areas,:) = [];
roiUsed = roiNames;
roiUsed(missing_areas) = [];

% add veno signal
X(size(X,1)+1,:) = (nanmean(avRespVeno_hr(t_select,:,:),3)');
roiUsed(length(roiUsed)+1) = {'Veno'};

% construct distance matrix
D_length = sqrt(sum(X.^2,2));
D_dot = (X*X');
D = D_dot./(D_length*D_length');
D(D>1) = 1;
D_mat = rad2deg(acos(D));

[Y,eigvals] = cmdscale(round(D_mat));
% figure,plot(eigvals)

% [Y] = mdscale(round(D_mat),2); % nonclassical

figure
plot(Y(:,1),Y(:,2),'.')
text(Y(:,1),Y(:,2),roiUsed,'FontSize',8)

figure('Position',[0 0 400 1000])
% imagesc(t_D,[1:length(roiUsed)],X,[-0.005 0.005])
imagesc(t_D,[1:length(roiUsed)],X,[-.1 .1])
set(gca,'YTick',[1:length(roiUsed)],'YTickLabel',roiUsed,'FontSize',9,'FontName','Arial')



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


%%
%% get distance matric for MDS, concatinate all subjects

t_select = find(t_hr<0.5);

% concatinate all subjects
X = reshape(avResp_hr(t_select,:,:),length(t_select),size(avResp_hr,2)*size(avResp_hr,3))';
t_D = t_hr(t_select);

% clear missing areas
missing_areas = isnan(X(:,1)) | sum(X,2)==0;
X(missing_areas,:) = [];
roiUsed = repmat(roiNames,1,size(avResp_hr,3));
roiUsed(missing_areas) = [];

% add veno signal
X = [X; squeeze(avRespVeno_hr(t_select,:,:))'];
roiUsed(length(roiUsed)+1:length(roiUsed)+size(avResp_hr,3)) = repmat({'Veno'},6,1);

% construct distance matrix
D_length = sqrt(sum(X.^2,2));
D_dot = (X*X');
D = D_dot./(D_length*D_length');
D(D>1) = 1;
D_mat = rad2deg(acos(D));

[Y,eigvals] = cmdscale(round(D_mat));
% figure,plot(eigvals)

% [Y] = mdscale(round(D_mat),2); % nonclassical
%%

% todo: add labels / remove subcortical?

% plot
figure,hold on
plot(Y(:,1),Y(:,2),'.')
text(Y(:,1),Y(:,2),roiUsed)

% add Veno colors:
these_areas = find(ismember(roiUsed,'Veno'));
plot(Y(these_areas,1),Y(these_areas,2),'r*','MarkerSize',20)

% add CSF colors:
these_areas = find(ismember(roiUsed,{'leftlateralventricle','rightlateralventricle'}));
plot(Y(these_areas,1),Y(these_areas,2),'g*','MarkerSize',20)

% add posterior Cerebral colors:
these_areas = find(ismember(roiUsed,{'rightparahippocampal','leftparahippocampal'}));
plot(Y(these_areas,1),Y(these_areas,2),'b*','MarkerSize',20)

% add anterior Cerebral colors:
these_areas = find(ismember(roiUsed,{'leftcaudalanteriorcingulate','rightcaudalanteriorcingulate'}));
plot(Y(these_areas,1),Y(these_areas,2),'m*','MarkerSize',20)

% add middle Cerebral colors:
these_areas = find(ismember(roiUsed,{'rightinsula','leftinsula'}));
plot(Y(these_areas,1),Y(these_areas,2),'c*','MarkerSize',20)

% add far away from main branches:
these_areas = find(ismember(roiUsed,{'rightinferiortemporal','leftinferiortemporal'}));
plot(Y(these_areas,1),Y(these_areas,2),'y*','MarkerSize',20)

% % add CSF:
% these_areas = find(ismember(roiUsed,'CSF'));
% plot(Y(these_areas,1),Y(these_areas,2),'y*','MarkerSize',20)

% % add ITG:
% these_areas = find(ismember(roiUsed,'rightinferiortemporal'));
% plot(Y(these_areas,1),Y(these_areas,2),'y*','MarkerSize',20)

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir '/figures/segmentation/allSubs_MDStest'])
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/allSubs_MDStest'])

%%
figure('Position',[0 0 400 1000])
imagesc(t_D,[1:length(roiUsed)],X,[-0.005 0.005])
set(gca,'YTick',[1:length(roiUsed)],'YTickLabel',roiUsed,'FontSize',9,'FontName','Arial')

%% plot some responses close and further from main branch of arteries

figure('Position',[0 0 300 800])

subplot(5,1,1),hold on
plot(t_hr,zeros(size(t_hr)),'Color',[.7 .7 .7])
plot([0 0],[-0.2 0.2],'Color',[.7 .7 .7])
title('posterior cerebral branch')
this_area = 45; % right parahippocampal
plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'k','LineWidth',1)

% this_area = 42; % right lingual
% plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'Color',[.5 .5 .5])

this_area = 38; % right ITG
plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'Color',[.5 .5 1])
xlim([t_hr(1) t_hr(end)])

subplot(5,1,2),hold on
plot(t_hr,zeros(size(t_hr)),'Color',[.7 .7 .7])
plot([0 0],[-0.2 0.2],'Color',[.7 .7 .7])

title('anterior cerebral branch')
this_area = 1; % caudal anterior cingulate 1 
plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'k','LineWidth',1)

% this_area = 52; % posterior cingulate 21 
% plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'Color',[.5 .5 .5])

this_area = 39; % isthmus cingulate 8 39
plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'Color',[.5 .5 1])

subplot(5,1,3),hold on
plot(t_hr,zeros(size(t_hr)),'Color',[.7 .7 .7])
plot([0 0],[-0.2 0.2],'Color',[.7 .7 .7])

title('middle cerebral branch')
this_area = 62; % right insula
plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'k','LineWidth',1)

this_area = 44; % middle temporal
plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'Color',[.5 .5 1])

% this_area = 56; % right rostral middle frontal
% plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'g')

subplot(5,1,4),hold on
plot(t_hr,zeros(size(t_hr)),'Color',[.7 .7 .7])
plot([0 0],[-0.2 0.2],'Color',[.7 .7 .7])

title('Veins')
plot(t_hr,squeeze(avRespVeno_hr(:,1,:)),'k','LineWidth',1)


subplot(5,1,5),hold on
title('Lateral ventricles')
plot(t_hr,zeros(size(t_hr)),'Color',[.7 .7 .7])
plot([0 0],[-0.2 0.2],'Color',[.7 .7 .7])
% average left and lateral ventricles 74 89
plot(t_hr,(squeeze(avResp_hr(:,74,:))+squeeze(avResp_hr(:,89,:)))/2,'k','LineWidth',1)
% Caudate:
% plot(t_hr,squeeze(avResp_hr(:,79,:)),'Color',[1 .5 .5])
% plot(t_hr,squeeze(avResp_hr(:,94,:)),'Color',[.5 1 .5])
% Putamen:
% plot(t_hr,squeeze(avResp_hr(:,80,:)),'Color',[1 .5 .5])
% plot(t_hr,squeeze(avResp_hr(:,95,:)),'Color',[.5 1 .5])

% % right lateral ventricle
% this_area = 89;
% plot(t_hr,squeeze(avResp_hr(:,this_area,:)),'r')

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir '/figures/segmentation/allSubs_tracessegm01'])
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/allSubs_tracessegm01'])

