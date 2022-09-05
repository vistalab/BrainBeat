clear all
close all

% script produced transformation matrix between the anatomy and functional
% scans, walk through step-by-step

[~,dDir] = bbPath;

%% Coregistration functionals to T1
 
sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,[1 2],[1 2],[1 2],[1 2]};


for ss = 1 % subject number loop
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};

    % Get the anatomicals:
    t1_name = fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'anat',...
        ['sub-' sub_label '_ses-' ses_label '_T1w.nii']);
    niAnatomy = niftiRead(t1_name);

    %%%%% coregister the functionals to the T1:
    for rr = 1:length(run_nrs{ss}) % runs

        run_nr = run_nrs{ss}(rr);
        
        % read functional data
        fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
            ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
        fmri_name = fullfile(dDir,fmri_BIDSname);
        ni1 = niftiRead(fmri_name);
        
        % directory and base name of how to save the coregitsration matrix
        save_dir = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label]);
        save_name_base = ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)];
        
        % select certain scans to align in ni1 (first, average)
        ni1 = ni;
        %%%%% use the first nifti to align, this one has the most structural info:
        % this worked best for subject 1, scan 2
        ni1.data=ni1.data(:,:,:,1); % use this
        %%%%% use the average functionals to align
%         ni1.data = mean(ni1.data(:,:,:,5:end),4); 
        
        ni1.dim = ni1.dim(1:3);
        ni1.pixdim = ni1.pixdim(1:3);
        
        niAnatomy.pixdim = niAnatomy.pixdim(1:3); % only uses the first three dimensions

        % align functionals to the T1:
        acpcXform = dtiRawAlignToT1(ni1,niAnatomy,[], [], false, 1); % last 1 adds a figure
        % this saves the realignment matrix in the folder of the functionals, 
        
        % Only if the current acpcXform is good enough, safe for use
        acpcXform_new = acpcXform;
        save(fullfile(save_dir,[save_name_base '_AcpcXform_new.mat']),'acpcXform_new')
    end
    
end


%% Check the coregistration

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,[1 2],[1 2],[1 2],[1 2]};

ss = 1; % subject number loop
sub_label = sub_labels{ss};
ses_label = ses_labels{ss};
acq_label = acq_labels{ss};
run_nr = run_nrs{ss}(1); % run 1

% Load T1
t1_name = fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'anat',...
    ['sub-' sub_label '_ses-' ses_label '_T1w.nii']);
niAnatomy = niftiRead(t1_name);

% Load functionals
fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
    ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
fmri_name = fullfile(dDir,fmri_BIDSname);
ni = niftiRead(fmri_name);

% Load coregistration
save_dir = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label]);
save_name_base = ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)];
load(fullfile(save_dir,[save_name_base '_AcpcXform_new.mat']),'acpcXform_new')


% Settings for overlay:
curPos = [1,18,35]; 
sliceThisDim = 1; %, z, y, x
imDims = [-90 -120 -120; 90 130 120];

%%%% Overlay 1: functionals and anatomy
niFunc = ni;
niFunc.data = mean(ni.data(:,:,:,1),4); % overlay the first functional
% niFunc.data = mean(ni.data(:,:,:,5:end),4); % overlay the mean functional
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform_new,sliceThisDim,imDims,curPos);



%%
%% %% Overlay anatomical/functionals and timeseries for fun
%%

%% The T2* data are here.  

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {1,1,[1 2],[1 2],[1 2],[1 2]};

ss = 1; % subject number loop
sub_label = sub_labels{ss};
ses_label = ses_labels{ss};
acq_label = acq_labels{ss};
run_nr = run_nrs{ss}(1); % run 1

% Load T1
t1_name = fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'anat',...
    ['sub-' sub_label '_ses-' ses_label '_T1w.nii']);
niAnatomy = niftiRead(t1_name);

% Load functionals
fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
    ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
fmri_name = fullfile(dDir,fmri_BIDSname);
ni = niftiRead(fmri_name);

% Load coregistration
save_dir = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label]);
save_name_base = ['sub-' sub_label '_ses-' ses_label '_task-rest_acq-' acq_label '_run-' int2str(run_nr)];
load(fullfile(save_dir,[save_name_base '_AcpcXform_new.mat']),'acpcXform_new')


% load time series and associated time
load(fullfile(save_dir,[save_name_base '_PPGtrigResponseT']),'t')
ppgTSname = fullfile(save_dir,[save_name_base '_PPGtrigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series


% load correlation between even and odd scans for colors of timeseries
ppgRname = fullfile(save_dir,[save_name_base '_codPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

%%%% Overlay 2: timeseries and anatomy
ppgTSplot = ppgTS;
ppgTSplot.data(:,:,:,t<-.1 | t>1)=[]; % plot these times from curve

% Scale time series amplitude by R, plots are generated with respect to the maximum.
niColor = ppgR; % use R2 map to scale later: r2_scale=sqrt(imgSlice2.^2);
maxTS = max(abs(ppgTSplot.data),[],4); % get the max of each curve
ppgTSplot.data = bsxfun(@rdivide,ppgTSplot.data,maxTS); % devide by the max of each curve (sets all curves to 1 max)
ppgTSplot.data = bsxfun(@times,ppgTSplot.data,niColor.data); % multiply by r^2 to set less reliable curves to zero

bbOverlayTimeseriesAnat(ppgTSplot,niColor,niAnatomy,acpcXform_new,sliceThisDim,imDims,curPos)

