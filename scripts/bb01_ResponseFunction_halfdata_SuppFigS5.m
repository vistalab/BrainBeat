
%% How well do we do with less data?

% this script runs through code to test whether less data would give
% similar estimated curves


% dependency: spm (spm_bsplinc)


%% Set base data directory 

clear all
close all

[~,dDir] = bbPath;


%% Preprocess increasingly more data and get PPG time series and calculate COD

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

% We decrease the number of scans in the subset :
scans_in_subsets = [879 720 480 240 120]; % all data to less

for ss = 1:length(sub_labels) % subjects/ses/acq
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    
    for rr = 1:length(run_nrs{ss}) % runs
        
        run_nr = run_nrs{ss}(rr);
        
        fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
            ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
        fmri_name = fullfile(dDir,fmri_BIDSname);
        save_dir = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label]);
        if ~exist(save_dir,'dir')
            disp(['creating output directory ' save_dir])
            mkdir(save_dir)
        end
        save_name_base = ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)];

        if ~exist(fmri_name,'file')
            clear ni
            error('filename %s does not exist',fmri_name)
        end
        niAll = niftiRead(fmri_name);

        disp(['Calculating PPG triggered responses for ' fmri_BIDSname])        

        %%%% Run through subsets
        
        % Define output
        ppgR.data = zeros([size(ppgR_all.data) length(scans_in_subsets)]);
        for kk = 1:length(scans_in_subsets)
            disp(['subset ' int2str(kk)])
            % compute the PPG triggered response matrix for the whole brain
            % and save as a nifti: we always need to go back from to end,
            % since that is the assumption for the PPG
            ni = niAll;
            ni.data = niAll.data(:,:,:,[size(niAll.data,4)-scans_in_subsets(kk)+1:size(niAll.data,4)]);
            ni.dim = size(ni.data);
            [~,t,response_matrix_odd,response_matrix_even] = bbResponse2physio(ni);       
            % average of all odd heartbeats:
            ni_odd = ni; ni_odd.data = response_matrix_odd;
            % average of all even heartbeats:
            ni_even = ni; ni_even.data = response_matrix_even;

            % compute the COD between the even and odd repetitions of
            % the PPG triggered BOLD signal for the whole brain:
            [out_R_map] = bbCod2physio(ni_odd,ni_even);
            ppgR.data(:,:,:,kk) = out_R_map;
        end
        
        % save ppgR as nifti
        ni2 = niAll;
        ni2.data = ppgR.data;
        ni2.fname = [save_name_base '_codPPG_lessData.nii.gz'];
        ni2.descrip = [ni2.descrip ';scansInSets=' num2str(scans_in_subsets)];
        niftiWrite(ni2,fullfile(save_dir,ni2.fname))
        
    end
end

%% Start here to load data, when previous already run

scans_in_subsets = [879 720 480 240 120];

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};
out = [];
for ss = 1:length(sub_labels) % subjects/ses/acq
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    rr = 1;        
    run_nr = run_nrs{ss}(rr);
        
    % get the ppg COD from subsets
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);
    ppgR_subsets = niftiRead([save_name_base '_codPPG_lessData.nii.gz']);
    out(ss).data = ppgR_subsets.data;
end

%%
relVox = zeros(length(sub_labels),length(scans_in_subsets));
corr_HighLow = zeros(length(sub_labels),length(scans_in_subsets));
for ss = 1:length(sub_labels) % subjects/ses/acq
    for kk = 1:length(scans_in_subsets)
        relVox(ss,kk) = length(find(out(ss).data(:,:,:,kk)>.5))./length(find(out(ss).data(:,:,:,1)>.5));
        best_R = out(ss).data(:,:,:,1); best_R = best_R(:);
        this_R = out(ss).data(:,:,:,kk); this_R = this_R(:);
        corr_HighLow(ss,kk) = corr(best_R,this_R,'type','Spearman');
    end
end

%%
figure('Position',[0 0 220 300])
subplot(2,1,1),hold on
bar(scans_in_subsets*.25,mean(relVox,1),'FaceColor',[1 1 1])
errorbar(scans_in_subsets*.25,mean(relVox,1),std(relVox,1),'k.')
set(gca,'XTick',round(sort(scans_in_subsets*.25)))
ylabel('Relative #voxels with COD>50')

subplot(2,1,2),hold on
bar(scans_in_subsets*.25,mean(corr_HighLow,1),'FaceColor',[1 1 1])
errorbar(scans_in_subsets*.25,mean(corr_HighLow,1),std(corr_HighLow,1),'k.')
set(gca,'XTick',round(sort(scans_in_subsets*.25)))
ylabel('Spearman R')
xlabel('Duration (s)')

% set(gcf,'PaperPositionMode','auto')    
% print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group',['testLessData_n6']))
% print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['testLessData_n6']))


%% load data to slices with reliable voxels when using less data in subject 1

sub_labels = {'1'};%,'2','3','4','5','1'}; 
ses_labels = {'1'};%,'1','1','1','1','2'}; 
acq_labels = {'4mmFA48'};%,'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1]};%,[1],[1],[1],[1],[1]};

ss = 1;
sub_label = sub_labels{ss};
ses_label = ses_labels{ss};
acq_label = acq_labels{ss};
    
rr = 1;
run_nr = run_nrs{ss}(rr);

% Get the anatomicals:
t1w_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'anat',...
            ['sub-' sub_label '_ses-' ses_label '_T1w.nii.gz']);
niAnatomy = niftiRead(fullfile(dDir,t1w_BIDSname));

% get the ppg COD from subsets
save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
    ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);
ppgR_subsets = niftiRead([save_name_base '_codPPG_lessData.nii.gz']);

% load coregistration matrix (for the functionals):
load([save_name_base '_AcpcXform_new.mat']);
acpcXform = acpcXform_new; clear acpcXform_new

%% plot slices with reliable voxels when using less data in subject 1

for index_subset = 1:5

    curPos = [-10,1,-20]; 
    sliceThisDim = 1; 
    imDims=[-90 -120 -120; 90 130 90];
    overlayPlot = ppgR_subsets;
    overlayPlot.data = overlayPlot.data(:,:,:,index_subset);
    cod_th = 0.5;
    
    for kk = -2
        curPos(1) = kk;
        bbOverlayFuncAnat(overlayPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,1,1,cod_th);
        title(['slice ' int2str(kk) ', R^2>' num2str(cod_th,3) ' #scans ' int2str(scans_in_subsets(index_subset))])
        set(gcf,'PaperPositionMode','auto')    
        % print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['testLessData_sub-' int2str(ss) '_nscans-' int2str(scans_in_subsets(index_subset))]))
    end
end

%%
%% now analyze run 2 and compare to that to subset from run 1
%% only possible for subjects with 2 runs
%%

scans_in_subsets = [879 720 480 240 120];

sub_labels = {'3','4','5','1'}; 
ses_labels = {'1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1 2],[1 2],[1 2],[1 2]};
out = [];
for ss = 1:length(sub_labels) % subjects/ses/acq
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    rr = 1;        
    run_nr = run_nrs{ss}(rr);
        
    % get the ppg COD from subsets
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);
    acpx_xform1 = load([save_name_base '_AcpcXform_new.mat']);
    ppgR_subsets = niftiRead([save_name_base '_codPPG_lessData.nii.gz']);
    out(ss).data = ppgR_subsets.data;
    out(ss).xform1 = acpx_xform1.acpcXform_new;
    
    % get the ppg COD from second run
    rr = 2;        
    run_nr = run_nrs{ss}(rr);
    save_name_base = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label],...
        ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)]);
    acpx_xform2 = load([save_name_base '_AcpcXform_new.mat']);
    ppgR_2 = niftiRead([save_name_base '_codPPG.nii.gz']);
    out(ss).ppgR_2 = ppgR_2.data;
    out(ss).xform2 = acpx_xform2.acpcXform_new;
end

%% we need to get data in same space, we go to anatomy for both, but retain voxel size...

relVox = zeros(length(sub_labels),length(scans_in_subsets));
corr_HighLow = zeros(length(sub_labels),length(scans_in_subsets));

bb = [-90 -120 -60; 90 96 130];
mmPerVox = [4 4 4];
bSplineParams = [1 1 1 0 0 0];

for ss = 1:length(sub_labels) % subjects/ses/acq

    img1 = out(ss).data;
    xform1 = inv(out(ss).xform1);
    [newImg1] = mrAnatResliceSpm(img1, xform1, bb, mmPerVox, bSplineParams);

    img2 = out(ss).ppgR_2;
    xform2 = inv(out(ss).xform2);
    [newImg2] = mrAnatResliceSpm(img2, xform2, bb, mmPerVox, bSplineParams);

    for kk = 1:length(scans_in_subsets)
        relVox(ss,kk) = length(find(newImg1(:,:,:,kk)>.5))./length(find(newImg2(:,:,:)>.5));
        best_R = newImg2; best_R = best_R(:);
        this_R = newImg1(:,:,:,kk); this_R = this_R(:);
        excl_nan = isnan(best_R) | isnan(this_R);
        corr_HighLow(ss,kk) = corr(best_R(~excl_nan),this_R(~excl_nan),'type','Spearman');
    end
end

%%
figure('Position',[0 0 220 300])
subplot(2,1,1),hold on
% plot(scans_in_subsets*.25,relVox,'k.')
bar(scans_in_subsets*.25,mean(relVox,1),'FaceColor',[1 1 1])
errorbar(scans_in_subsets*.25,mean(relVox,1),std(relVox,1),'k.')
set(gca,'XTick',round(sort(scans_in_subsets*.25)))
ylabel('Relative #voxels with COD>50')

subplot(2,1,2),hold on
bar(scans_in_subsets*.25,mean(corr_HighLow,1),'FaceColor',[1 1 1])
errorbar(scans_in_subsets*.25,mean(corr_HighLow,1),std(corr_HighLow,1),'k.')
set(gca,'XTick',round(sort(scans_in_subsets*.25)))
ylim([0 1])
ylabel('Spearman R')
xlabel('Duration (s)')

% set(gcf,'PaperPositionMode','auto')    
% print('-painters','-r300','-depsc',fullfile(dDir,'derivatives','brainbeat','group',['testRetestData_n4']))
% print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group',['testRetestData_n4']))

