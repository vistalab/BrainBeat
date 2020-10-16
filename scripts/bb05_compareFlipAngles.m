clear all
close all

%% Script to plot different flip angles
%


%% The T2* data are here.  
% dDir = '/biac4/wandell/data/BrainBeat/data';
% dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

in_data = 'PPG';
s = 3;
s_info = bb_subs(s);
subj = s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

if s==1
    all_scans = [2 1 3];
elseif s==2
    all_scans = [2 1 3];
elseif s==3
    all_scans = [1 2 3];
end

im_Max = zeros(length(all_scans),1);

for kk = 1:length(all_scans)
    scan_nr = all_scans(kk);
    % Get the functionals:
    scan = s_info.scan{scan_nr};
    scanName = s_info.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
    ni = niftiRead(fmri);
    meanNi(kk).im = mean(ni.data(:,:,:,5:end),4);
    imMax(kk) = max(meanNi(kk).im(:));
end


figure('Position',[0 0 800 600])
for kk = 1:length(all_scans)
    subplot(2,length(all_scans),kk)
    imagesc(meanNi(kk).im(:,:,20)',[1 max(imMax)])
    axis image
    
    subplot(2,length(all_scans),kk+3)
    imagesc(squeeze(meanNi(kk).im(30,:,:))',[1 max(imMax)])
    axis image
    axis xy
%     clear ni meanNi1
end

%% Load a segmentation

in_data = 'PPG';
s = 3;
s_info = bb_subs(s);
subj = s_info.subj;

if s==1
    all_scans = [2 1 3];
elseif s==2
    all_scans = [2 1 3];
elseif s==3
    all_scans = [1 2 3];
%     all_scans = [4 5 6];
%     all_scans = [7 8 9];
end

figure('Position',[0 0 500 500]),hold on

for kk = 1:length(all_scans)

    scan_nr = all_scans(kk);

    % Load the functionals:
    scan = s_info.scan{scan_nr};
    scanName = s_info.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
    ni = niftiRead(fmri);
    meanNi1 = mean(ni.data(:,:,:,5:end),4);

    % Load the segmentation
    segName = fullfile(dDir,subj,scan,[scanName '_combineSegm.nii.gz']);
    seg = niftiRead(segName);

    % Load the correlation with heartbeat (made with bbCorrelate2physio):
    roiNames = {'GM','WM','Ventricles','CSF','Veins'};    
    roiColor = [.6 .6 .6;.95 .95 .95;0 1 0;0 .5 .5;0 0 1];
    ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
    ppgR = niftiRead(ppgRname); % correlation with PPG

    subplot(3,3,kk),hold on
    % try boxplot
%     roiColor = [.3 .3 .3;.8 .8 .8;0 1 0;0 .5 .5;0 0 1];
%     theseVals = [];
%     theseGroups = [];
%     for rr = 1:5
%         theseVals = [theseVals; meanNi1(seg.data==rr)];
%         theseGroups = [theseGroups rr*ones(1,length(find(seg.data==rr)))];
%     end
%     boxplot(theseVals,theseGroupsm'colors',roiColor); 
    
    for rr = 1:5
        bar(rr,mean(meanNi1(seg.data==rr)),'FaceColor',roiColor(rr,:)); 
        errorbar(rr,mean(meanNi1(seg.data==rr)),std(meanNi1(seg.data==rr)),'k')
    end
    set(gca,'XTick',[1:5],'XTickLabel',roiNames);
    xtickangle(45)
    title(['FA ' int2str(s_info.scanFA{scan_nr})])
    ylabel('mean signal')
    
    subplot(3,3,3+kk),hold on
    plot([0 6],[1000 1000;2000 2000]','k')
    for rr = 1:5
        bar(rr,sum(ppgR.data(seg.data==rr)>0.3),'FaceColor',roiColor(rr,:)); 
    end
    set(gca,'XTick',[1:5],'XTickLabel',roiNames);
    xtickangle(45)
    title(['FA ' int2str(s_info.scanFA{scan_nr})])
    ylabel('# voxels with R>0.3')
    ylim([0 3000])

    subplot(3,3,6+kk),hold on
    plot([0 6],[.25 .25;.5 .5;.75 .75]','k')
    for rr = 1:5
        bar(rr,length(find(ppgR.data(seg.data==rr)>0.3))./length(find(seg.data==rr)),'FaceColor',roiColor(rr,:)); 
    end
    set(gca,'XTick',[1:5],'XTickLabel',roiNames); 
    xtickangle(45)
    title(['FA ' int2str(s_info.scanFA{scan_nr})])
    ylabel('% voxels with R>0.3')
    ylim([0 1])
end

% 
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','flipAngle',['sub-' int2str(s) '_FAdifferences']))
% print('-painters','-r300','-depsc',fullfile(dDir,'figures','flipAngle',['sub-' int2str(s) '_FAdifferences']))
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','flipAngle',['sub-' int2str(s) '_FAdifferences_set2']))
% print('-painters','-r300','-depsc',fullfile(dDir,'figures','flipAngle',['sub-' int2str(s) '_FAdifferences_set2']))
% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,'figures','flipAngle',['sub-' int2str(s) '_FAdifferences_3mm']))
% print('-painters','-r300','-depsc',fullfile(dDir,'figures','flipAngle',['sub-' int2str(s) '_FAdifferences_3mm']))


%% plot T1 with segmentation overlay

in_data = 'PPG';
s = 2;
s_info = bb_subs(s);
subj = s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

% Load the functionals:
scan_nr = 3; % pick one functional, because segmentation is resliced into functional
scan = s_info.scan{scan_nr};
scanName = s_info.scanName{scan_nr};

% Load coregistration matrix:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

roiNames = {'GM','WM','Ventricles','CSF','Veins'};
roiColor = [.6 .6 .6;.95 .95 .95;0 1 0;0 .5 .5;0 0 1];

figure
sliceThisDim = 1;
if s==1
    imDims = [-90 -120 -100; 90 130 120];
    curPos = [4,18,5];
end
bbOverlayDotsAnat_PickColor(seg,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,roiColor);

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'figures','flipAngle',['sub-' int2str(s) '_scan-' int2str(scan_nr) '_FAdifferences']))


%% Plot as a function of FA to compare to prediction

in_data = 'PPG';
s = 1;
s_info = bb_subs(s);
subj = s_info.subj;

if s==1
    all_scans = [2 1 3];
elseif s==2
    all_scans = [2 1 3];
elseif s==3
%     all_scans = [1 2 3];
    all_scans = [4 5 6];
%     all_scans = [7 8 9];
end

thisFA = zeros(length(all_scans),1);
avVal = zeros(length(all_scans),5);

for kk = 1:length(all_scans)

    scan_nr = all_scans(kk);

    % Load the functionals:
    scan = s_info.scan{scan_nr};
    scanName = s_info.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
    ni = niftiRead(fmri);
    meanNi1 = mean(ni.data(:,:,:,5:end),4);

    % Load the segmentation
    segName = fullfile(dDir,subj,scan,[scanName '_combineSegm.nii.gz']);
    seg = niftiRead(segName);

    % Load the correlation with heartbeat (made with bbCorrelate2physio):
    roiNames = {'GM','WM','Ventricles','CSF','Veins'};    
    roiColor = [.6 .6 .6;.95 .95 .95;0 1 0;0 .5 .5;0 0 1];
    ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
    ppgR = niftiRead(ppgRname); % correlation with PPG

    thisFA(kk) = s_info.scanFA{scan_nr};
    
    for rr = 1:5
        avVal(kk,rr) = mean(meanNi1(seg.data==rr));
    end
end
%%
figure,hold on
roiColor = [.3 .3 .3;.7 .7 .7;0 1 0;0 .5 .5;0 0 1];
for rr = 1:5
    plot(thisFA,avVal(:,rr),'Color',roiColor(rr,:))
end