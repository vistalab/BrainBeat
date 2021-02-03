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
dkt_table_surface = readtable(fullfile(dDir,'dkt_areas_surface.tsv'),'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
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

    % resampling to hartbeat time with 128 timepoints t_hr
    t_temp = linspace(min(t_sel),max(t_sel),128);

    avResp_hr(:,:,ss) = interp1(t_sel,avResp_sel,t_temp);
    avRespVeno_hr(:,:,ss) = interp1(t_sel,avRespVeno_sel,t_temp);
    ppgCurve_hr(:,ss) = interp1(ppgt_sel,ppgCurve_sel,t_temp);
    
    clear t_sel ppgt_sel t_temp % t_sel and ppgt_sel are subject specific
    clear avResp_sel avRespVeno_sel ppgCurve_sel avRespVeno avResp ...
        ppgCurveT ppgCurve ppg_cycle ppg_onsets % all other subject specific vars
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

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir '/figures/segmentation/Fig1D_allSubs_tracessegm01'])
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/Fig1D_allSubs_tracessegm01'])

%% check render

% Select a subject number and hemisphere
ss = 1;
hemi = 'l'; % 'l' or 'r'

s_nr = all_subs(ss); %[1 2 3 4 5 6]
subs = bb_subs(s_nr);
subj = subs.subj;
% load gifti
gi = gifti(fullfile(dDir,subj,subs.anat,['T1w_' hemi 'h_render.gii']));

% load DKT surface labels
surface_labels_name = fullfile(dDir,subj,'freesurfer','label',[hemi 'h.aparc.DKTatlas.annot']);

% surface_labels = MRIread(surface_labels_name);
[vertices, label, colortable] = read_annotation(surface_labels_name);

% mapping labels to colortable, because labels are strange (e.g. 9182740)
vert_label = label; 
for kk = 1:size(colortable.table,1) % 36 labels
    vert_label(label==colortable.table(kk,5)) = kk;
end
clear label 

% vert_label now corresponds to dkt_table_surface

% make a colormap for the labels
cmap = [dkt_table_surface.r1 dkt_table_surface.g1 dkt_table_surface.b1];
% replace NaNs in colormap with gray
cmap(isnan(cmap)) = 0.7;

figure
ieeg_RenderGiftiLabels(gi,vert_label,cmap,colortable.struct_names)
ieeg_viewLight(90,0)
%%
% create colors for rendering and plots:
nr_vals = 4;
cm1 = [[0.7:0.3/(nr_vals-1):1]' [0.7:-0.2/(nr_vals-1):0.5]' [0.7:-0.2/(nr_vals-1):0.5]'];
figure,imagesc(1:100)
colormap(cm1)
nr_vals = 4;
cm2 = [[0.7:-0.2/(nr_vals-1):0.5]' [0.7:0.3/(nr_vals-1):1]' [0.7:-0.2/(nr_vals-1):0.5]'];
nr_vals = 5;
cm3 = [[0.7:-0.2/(nr_vals-1):0.5]' [0.7:-0.2/(nr_vals-1):0.5]' [0.7:0.3/(nr_vals-1):1]'];

%% plot responses in these colors

% we have avResp_hr with 106 regions
% there are 31 surface areas with dkt_table.DKT_nr corresponding to dkt_table_surface.DKT_nr
% so we use the colors from dkt_table_surface.DKT_nr

figure
subplot(3,1,1),hold on
ind_art = 1001:1004; % main branch of anterior cerebral artery
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    plot(t_hr,mean(avResp_hr(:,these_responses,:),3),'Color',these_colors,'LineWidth',2)
end
ylim([-0.3 0.3])

subplot(3,1,2),hold on
ind_art = 2001:2004; % main branch of middle cerebral artery
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    plot(t_hr,mean(avResp_hr(:,these_responses,:),3),'Color',these_colors,'LineWidth',2)
end
ylim([-0.3 0.3])

subplot(3,1,3),hold on
ind_art = 3001:3005; % main branch of posterior cerebral artery
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    plot(t_hr,mean(avResp_hr(:,these_responses,:),3),'Color',these_colors,'LineWidth',2)
end
ylim([-0.3 0.3])
xlabel('heartbeat cycle')