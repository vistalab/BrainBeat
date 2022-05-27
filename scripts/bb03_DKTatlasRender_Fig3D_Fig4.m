
% This script generates Figures 3D and 4 from the manuscript titled:
%
% Measuring brain beats: cardiac-aligned fast fMRI signals
% Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell
%

clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%% Get the average responses across subjects

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,1,1,1,1};

dkt_table = readtable('dkt_areas.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dkt_table_surface = readtable('dkt_areas_surface.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
roiNames = dkt_table.label;
roiCodes = dkt_table.label_nr;

t_hr = linspace(-.5,1.5,128);
avResp_hr = NaN(length(t_hr),length(roiNames),length(sub_labels));
avRespVeno_hr = NaN(length(t_hr),1,length(sub_labels));
ppgCurve_hr = NaN(length(t_hr),length(sub_labels));
    
for ss = 1:length(sub_labels)
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};

    rr = 1;% just use first run
    run_nr = run_nrs{ss}(rr);
    
    % Functionals
    fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
    fmri_name = fullfile(dDir,fmri_BIDSname);
    ni = niftiRead(fmri_name);

    % get physio stuff we need:
    physio      = physioCreate('nifti',ni);
    ppgCurve    = physioGet(physio,'ppg ppgcurve');
    ppgCurveT   = physioGet(physio,'ppg ppgtcurve');
    srate       = 1/bbGet(ni,'tr');

    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);

    ppgResp = niftiRead([save_name_base '_PPGtrigResponse.nii.gz']); % ppg triggered time series
    ppgT = load([save_name_base '_PPGtrigResponseT.mat'],'t');

    %%%% Scale the time-series matrix by the reliability
    % Get odd/even COD (made with bbCod/Correlate2physio):
    ppgR = niftiRead([save_name_base '_codPPG.nii.gz']); % COD between even and odd heartbeats
    
    % load segmentation
    niSegm = niftiRead([save_name_base '_r_DKTatlas_aseg.nii.gz']);
    % matrix to vector
    segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);

    % also get the Venogram segmentation from niSegm2
    % niSegm2 has labels 1:5 with names {'GM','WM','Ventricles','CSF','Veno'};
    niSegm2 = niftiRead([save_name_base '_combineSegm.nii.gz']);
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
    ppg_cycle = 1./physioGet(physio,'PPGrate');
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
%% Figure 3D: plot set of responses

n_subs = size(avResp_hr,3);

figure('Position',[0 0 150 280])

subplot(4,1,1),hold on
this_area = 1; % caudal anterior cingulate 1 
% title(roiNames{this_area});
up_ci = mean(avResp_hr(:,this_area,:),3) + 2*std(avResp_hr(:,this_area,:),[],3)/sqrt(n_subs);
low_ci = mean(avResp_hr(:,this_area,:),3) - 2*std(avResp_hr(:,this_area,:),[],3)/sqrt(n_subs);
fill([t_hr t_hr(end:-1:1)],100*[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_hr,100*squeeze(avResp_hr(:,this_area,:)),'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,1,2),hold on
% title('Lateral ventricles')
thisSignal = (squeeze(avResp_hr(:,74,:))+squeeze(avResp_hr(:,89,:)))/2;
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
fill([t_hr t_hr(end:-1:1)],100*[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% average left and lateral ventricles 74 89
plot(t_hr,100*thisSignal,'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

subplot(4,1,3),hold on
% title('Veins')
up_ci = mean(avRespVeno_hr(:,1,:),3) + 2*std(avRespVeno_hr(:,1,:),[],3)/sqrt(n_subs);
low_ci = mean(avRespVeno_hr(:,1,:),3) - 2*std(avRespVeno_hr(:,1,:),[],3)/sqrt(n_subs);
fill([t_hr t_hr(end:-1:1)],100*[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_hr,100*squeeze(avRespVeno_hr(:,1,:)),'Color',[.5 .5 .5],'LineWidth',1)
xlim([t_hr(1) t_hr(end)])

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

set(gcf,'PaperPositionMode','auto')    
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
    ['Figure3D_allsubs']))
print('-r300','-depsc2',fullfile(dDir,'derivatives','figures',...
    ['Figure3D_allsubs']))

%% Figure 4A-B: Render subject 1 with DKT colored by arterial branch

% Select a subject number and hemisphere
ss = 1;
hemi = 'l'; % 'l' or 'r'

% Functionals
sub_label = sub_labels{ss};
ses_label = ses_labels{ss};
acq_label = acq_labels{ss};

% load gifti
gi = gifti(fullfile(dDir,'derivatives','surfaces',['sub-' sub_label],['ses-' ses_label],...
    ['T1w_' hemi 'h_render.gii']));

% load DKT surface labels
surface_labels_name = fullfile(dDir,'derivatives','freesurfer',['sub-' sub_label],['ses-' ses_label],...
    'label',[hemi 'h.aparc.DKTatlas.annot']);

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

set(gcf,'PaperPositionMode','auto')    
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
    ['Figure4A_sub-' sub_label '_renderArterialZones_mesial']))
ieeg_viewLight(270,0)
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
    ['Figure4B_sub-' sub_label '_renderArterialZones_lateral']))


%% Figure 4A-B colormap
cm1 = [[dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==1001,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==1001,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==1001,1))];...
    [dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==1002,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==1002,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==1002,1))];...
    [dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==1003,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==1003,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==1003,1))];...
    [dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==1004,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==1004,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==1004,1))]];
cm2 = [[dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==2001,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==2001,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==2001,1))];...
    [dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==2002,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==2002,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==2002,1))];...
    [dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==2003,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==2003,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==2003,1))];...
    [dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==2004,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==2004,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==2004,1))]];
cm3 = [[dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==3001,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==3001,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==3001,1))];...
    [dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==3002,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==3002,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==3002,1))];...
    [dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==3003,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==3003,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==3003,1))];...
    [dkt_table_surface.r1(find(dkt_table_surface.ind_arterial==3004,1)) dkt_table_surface.g1(find(dkt_table_surface.ind_arterial==3004,1)) dkt_table_surface.b1(find(dkt_table_surface.ind_arterial==3004,1))]];

figure('Position',[0 0 200 50]),
imagesc(1:100)
colormap(cm1)
axis off
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
    ['Figure4A_cm1']))

figure('Position',[0 0 200 50]),
imagesc(1:100)
colormap(cm2)
axis off
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
    ['Figure4A_cm2']))

figure('Position',[0 0 200 50]),
imagesc(1:100)
colormap(cm3)
axis off
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
    ['Figure4A_cm3']))


%% Figure 4C: plot responses in these colors

% we have avResp_hr with 106 regions
% there are 31 surface areas with dkt_table.DKT_nr corresponding to dkt_table_surface.DKT_nr
% so we use the colors from dkt_table_surface.DKT_nr

figure('Position',[0 0 200 400])
subplot(3,1,1),hold on
ind_art = 1001:1004; % main branch of anterior cerebral artery
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    plot(t_hr,100*mean(avResp_hr(:,these_responses,:),3),'Color',these_colors,'LineWidth',2)
end
ylim([-1 1])

subplot(3,1,2),hold on
ind_art = 2001:2004; % main branch of middle cerebral artery
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    plot(t_hr,100*mean(avResp_hr(:,these_responses,:),3),'Color',these_colors,'LineWidth',2)
end
ylim([-1 1])

subplot(3,1,3),hold on
ind_art = 3001:3004; % main branch of posterior cerebral artery
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    plot(t_hr,100*mean(avResp_hr(:,these_responses,:),3),'Color',these_colors,'LineWidth',2)
end
ylim([-1 1])
xlabel('heartbeat cycle')

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
    ['Figure4C_Artery_Waveforms']))
print('-r300','-depsc',fullfile(dDir,'derivatives','figures',...
    ['Figure4C_Artery_Waveforms']))

%% Figure 4D: get amplitude across subjects and test for significance

t_int = t_hr>-0.5 & t_hr<0.5; % time for max-min contrast

figure('Position',[0 0 200 400])
subplot(3,1,1),hold on
ind_art = 1001:1004; % main branch of anterior cerebral artery
art1_amp = zeros(length(ind_art),6,length(t_hr)); % 6 subjects
these_colors1 = zeros(length(ind_art),3);
for kk = 1:length(ind_art) 
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors1(kk,:) = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    art1_amp(kk,:,:) = 100*squeeze(mean(avResp_hr(:,these_responses,:),2))'; % average across multiple areas 
end

for kk = 1:length(ind_art)     
    % take contrast across a section
    bar(kk,max(mean(art1_amp(kk,:,t_int)),[],3) - min(mean(art1_amp(kk,:,t_int)),[],3),'FaceColor',these_colors1(kk,:))
    plot(kk,max(art1_amp(kk,:,t_int),[],3) - min(art1_amp(kk,:,t_int),[],3),'k.')

end
plot(max(art1_amp(:,:,t_int),[],3) - min(art1_amp(:,:,t_int),[],3),'LineWidth',1)
set(gca,'XTick',[1:4])

subplot(3,1,2),hold on
ind_art = 2001:2004; % main branch of middle cerebral artery
art2_amp = zeros(length(ind_art),6,length(t_hr)); % 6 subjects
these_colors2 = zeros(length(ind_art),3);
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors2(kk,:) = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    art2_amp(kk,:,:) = 100*squeeze(mean(avResp_hr(:,these_responses,:),2))'; % average across multiple areas 
end

for kk = 1:length(ind_art) 
    % take contrast across a section
    bar(kk,max(mean(art2_amp(kk,:,t_int)),[],3) - min(mean(art2_amp(kk,:,t_int)),[],3),'FaceColor',these_colors2(kk,:))
    plot(kk,max(art2_amp(kk,:,t_int),[],3) - min(art2_amp(kk,:,t_int),[],3),'k.')
end
plot(max(art2_amp(:,:,t_int),[],3) - min(art2_amp(:,:,t_int),[],3),'LineWidth',1)
set(gca,'XTick',[1:4])

subplot(3,1,3),hold on
ind_art = 3001:3004; % main branch of posterior cerebral artery
art3_amp = zeros(length(ind_art),6,length(t_hr)); % 6 subjects
these_colors3 = zeros(length(ind_art),3);
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors3(kk,:) = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    art3_amp(kk,:,:) = 100*squeeze(mean(avResp_hr(:,these_responses,:),2))'; % average across multiple areas 
end

for kk = 1:length(ind_art) 
    % take contrast across a section
    bar(kk,max(mean(art3_amp(kk,:,t_int)),[],3) - min(mean(art3_amp(kk,:,t_int)),[],3),'FaceColor',these_colors3(kk,:))
    plot(kk,max(art3_amp(kk,:,t_int),[],3) - min(art3_amp(kk,:,t_int),[],3),'k.')
end
plot(max(art3_amp(:,:,t_int),[],3) - min(art3_amp(:,:,t_int),[],3),'LineWidth',1)
set(gca,'XTick',[1:4])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','figures',...
    ['Figure4D_Artery_Contrast']))
print('-r300','-depsc',fullfile(dDir,'derivatives','figures',...
    ['Figure4D_Artery_Contrast']))

