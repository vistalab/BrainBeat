
% This script makes Figure 8 from:
%
%%%%% Measuring brain beats: cardiac-aligned fast fMRI signals %%%%%
% Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell
%
%

clear all
close all

%% Data directory 

[~,dDir] = bbPath;

%% Get the average responses across subjects

sub_labels = {'09','10','21','29','32'}; 

dkt_table = readtable('dkt_areas.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
dkt_table_surface = readtable('dkt_areas_surface.tsv','FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
roiNames = dkt_table.label;
roiCodes = dkt_table.label_nr;

avResp_all = NaN(248,length(roiNames),length(sub_labels));
% ppgCurve_hr = NaN(248,length(sub_labels));
%     
for kk = 1:length(sub_labels)
     
    sub_label = ['sub-IRIS0' sub_labels{kk}];
    
    % load cardiac gated timeseries
    ppgResp = niftiRead(fullfile(dDir,'derivatives','rs_fmri',sub_label,...
        ['avgGatedSignalNormalizedPercentage.nii.gz']));

    % load times
    t = load(fullfile(dDir,'derivatives','rs_fmri',sub_label,...
        ['avgTimes.mat']));
    
    % load segmentation
    niSegm = niftiRead(fullfile(dDir,'derivatives','rs_fmri',sub_label,...
        ['fs_aparc_resampled.nii.gz']));
    % matrix to vector
    segmVect = reshape(niSegm.data,[size(niSegm.data,1) * size(niSegm.data,2) * size(niSegm.data,3)],1);
    
    % put the data in a matrix Voxel X Time
    respMat = reshape(ppgResp.data,[size(ppgResp.data,1) * size(ppgResp.data,2) * size(ppgResp.data,3)],size(ppgResp.data,4));

    for ll = 1:length(roiNames)
        avResp_all(:,ll,kk) = median(respMat(ismember(segmVect,roiCodes(ll)),:),1);
    end
    
end

%%
%% Figure 4D: plot set of responses
t_hr = t.interpTimesNormalizedC;

% devide by 1000 and subtract a second to match SMS plots
t_match = t_hr/1000-1;

n_subs = size(avResp_all,3);

figure('Position',[0 0 150 280])

subplot(2,1,1),hold on
this_area = 62; % insula
% title(roiNames{this_area});
up_ci = mean(avResp_all(:,this_area,:),3) + 2*std(avResp_all(:,this_area,:),[],3)/sqrt(n_subs);
low_ci = mean(avResp_all(:,this_area,:),3) - 2*std(avResp_all(:,this_area,:),[],3)/sqrt(n_subs);
fill([t_match t_match(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
plot(t_match,squeeze(avResp_all(:,this_area,:)),'Color',[.5 .5 .5],'LineWidth',1)
xlim([-0.5 1.5])

subplot(2,1,2),hold on
% title('Lateral ventricles')
thisSignal = (squeeze(avResp_all(:,74,:))+squeeze(avResp_all(:,89,:)))/2;
up_ci = mean(thisSignal,2) + 2*std(thisSignal,[],2)/sqrt(n_subs);
low_ci = mean(thisSignal,2) - 2*std(thisSignal,[],2)/sqrt(n_subs);
fill([t_match t_match(end:-1:1)],[up_ci; low_ci(end:-1:1)],[0 0 0],'EdgeColor',[0 0 0])
% average left and lateral ventricles 74 89
plot(t_match,thisSignal,'Color',[.5 .5 .5],'LineWidth',1)
xlim([-0.5 1.5])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-depsc',[dDir '/derivatives/figures/Fig9AB_allSubs_slowfMRIventricleInsula'])
print('-painters','-r300','-dpng',[dDir '/derivatives/figures/Fig9AB_allSubs_slowfMRIventricleInsula'])


%%
% create colors for rendering and plots:
% nr_vals = 4;
% cm1 = [[.7 .7 .7]; 0.85*ones(nr_vals,1) [0.8:-0.8/(nr_vals-1):0]' [0.8:-0.8/(nr_vals-1):0]'];
% nr_vals = 4;
% cm2 = [[.7 .7 .7]; [0.8:-0.8/(nr_vals-1):0]' 0.85*ones(nr_vals,1) [0.8:-0.8/(nr_vals-1):0]'];
% nr_vals = 4;
% cm3 = [[.7 .7 .7]; [0.8:-0.8/(nr_vals-1):0]' [0.8:-0.8/(nr_vals-1):0]' 0.85*ones(nr_vals,1)];
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
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/Fig3A_s' int2str(s_nr) '_cm1'])

figure('Position',[0 0 200 50]),
imagesc(1:100)
colormap(cm2)
axis off
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/Fig3A_s' int2str(s_nr) '_cm2'])

figure('Position',[0 0 200 50]),
imagesc(1:100)
colormap(cm3)
axis off
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-depsc',[dDir '/figures/segmentation/Fig3A_s' int2str(s_nr) '_cm3'])


%% Figure 5B: plot responses in these colors

% translate t_hr to timing matching fast SMS sequence
% devide by 1000 and subtract a second
t_match = t_hr/1000-1;

% we have avResp_all with 106 regions
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
    plot(t_match,mean(avResp_all(:,these_responses,:),3),'Color',these_colors,'LineWidth',2)
end
xlim([-.5 1.5])

subplot(3,1,2),hold on
ind_art = 2001:2004; % main branch of middle cerebral artery
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    plot(t_match,mean(avResp_all(:,these_responses,:),3),'Color',these_colors,'LineWidth',2)
end
xlim([-.5 1.5])

subplot(3,1,3),hold on
ind_art = 3001:3004; % main branch of posterior cerebral artery
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    plot(t_match,mean(avResp_all(:,these_responses,:),3),'Color',these_colors,'LineWidth',2)
end
xlim([-.5 1.5])
xlabel('heartbeat cycle')

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-depsc',[dDir '/derivatives/figures/Fig9C_ArteryAvgs_slowfMRI'])
print('-painters','-r300','-dpng',[dDir '/derivatives/figures/Fig9C_ArteryAvgs_slowfMRI'])

%% Figure 5C: get amplitude across subjects and test for significance

% t_int = t_hr>-0.2 & t_hr<0; % for mean
t_int = t_hr>0 & t_hr<1000; % for max-min contrast

figure('Position',[0 0 200 400])
subplot(3,1,1),hold on
ind_art = 1001:1004; % main branch of anterior cerebral artery
art1_amp = zeros(length(ind_art),5,length(t_hr)); % 6 subjects
these_colors1 = zeros(length(ind_art),3);
for kk = 1:length(ind_art) 
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors1(kk,:) = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    art1_amp(kk,:,:) = squeeze(mean(avResp_all(:,these_responses,:),2))'; % average across multiple areas 
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
art2_amp = zeros(length(ind_art),5,length(t_hr)); % 5 subjects
these_colors2 = zeros(length(ind_art),3);
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors2(kk,:) = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    art2_amp(kk,:,:) = squeeze(mean(avResp_all(:,these_responses,:),2))'; % average across multiple areas 
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
art3_amp = zeros(length(ind_art),5,length(t_hr)); % 5 subjects
these_colors3 = zeros(length(ind_art),3);
for kk = 1:length(ind_art)
    % find areas with first ind_art
    dkt_surface_index = find(dkt_table_surface.ind_arterial == ind_art(kk));
    these_dkt = dkt_table_surface.DKT_nr(dkt_surface_index);
    these_colors3(kk,:) = [dkt_table_surface.r1(dkt_surface_index(1)) dkt_table_surface.g1(dkt_surface_index(1)) dkt_table_surface.b1(dkt_surface_index(1))];
    % find matching volume index from dkt_table
    these_responses = ismember(dkt_table.DKT_nr,these_dkt);
    art3_amp(kk,:,:) = squeeze(mean(avResp_all(:,these_responses,:),2))'; % average across multiple areas 
end

for kk = 1:length(ind_art) 
    % take contrast across a section
    bar(kk,max(mean(art3_amp(kk,:,t_int)),[],3) - min(mean(art3_amp(kk,:,t_int)),[],3),'FaceColor',these_colors3(kk,:))
    plot(kk,max(art3_amp(kk,:,t_int),[],3) - min(art3_amp(kk,:,t_int),[],3),'k.')
end
plot(max(art3_amp(:,:,t_int),[],3) - min(art3_amp(:,:,t_int),[],3),'LineWidth',1)
set(gca,'XTick',[1:4])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-depsc',[dDir '/derivatives/figures/Fig9D_ArteryAvgs_bar_slowfMRI'])
print('-painters','-r300','-dpng',[dDir '/derivatives/figures/Fig9D_ArteryAvgs_bar_slowfMRI'])

