clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses
%
%

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
dDir = '/biac4/wandell/data/BrainBeat/data';

% chdir(dDir)

%% Subject info

subj_nr = 2;
s_info = bb_subs(subj_nr);
subj=s_info.subj;


%% Get the anatomicals:

niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

%% Get the MRVenogram:

niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii.gz']));

% get Xform  venogram to the T1:
xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));

%%
%% plot 3 flip angles at the same time
%%

if subj_nr == 2
    scans_use=[2 1 3];
    fa_nr = [25 34 48]; 
    max_plot = 300;
elseif subj_nr == 3
    scans_use=[2 1 3];
    fa_nr = [20 36 48]; 
    max_plot = 100;
end

sl_plot = 20;
figure('Position',[0 0 800 300])

for s=1:3 % scan number
    subplot(1,3,s)

    % load functional data
    scan_nr = scans_use(s);
    scan=s_info.scan{scan_nr};
    scanName=s_info.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
    ni = niftiRead(fmri);

    % % load coregistration matrix:
    % load(fullfile(dDir,subj,scan,[scanName 'AcpcXform.mat']))
    
    % create mean of scans 5:end
    mean_nii = mean(ni.data(:,:,:,5:end),4);
    
    imagesc(mean_nii(:,:,sl_plot)',[0 max_plot])
    axis image
    title(['fa = ' int2str(fa_nr(s)) ' max plot = ' int2str(max_plot)])

end

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',['./figures/check_fa/' s_info.subj ])
