clear all
close all

% script produced transformation matrix between the anatomy and functional
% scans, walk through step-by-step

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
% dDir = '/biac4/wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% chdir(dDir)

%% coregistration functionals to T1 (course)
% this can work well for FA=25, but not for FA=36/48, thus: allign FA=25,
% and align 36/48 to FA=25 with code further below
 
for s = 1
    s_info = bb_subs(s);
    subj=s_info.subj;

    % Get the anatomicals:
    niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

%     % Get the MRVenogram:
%     niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));

    %%%%% coregister the functionals to the T1:

    for scan_nr = 1
        scan=s_info.scan{scan_nr};
        scanName=s_info.scanName{scan_nr};

        fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
        ni = niftiRead(fmri);

        ni1=ni;
        %%%%% use the first nifti to allign, this one has the most structural info:
%         ni1.data=ni1.data(:,:,:,1); % use this
        %%%%% use the average functionals to allign
        ni1.data=mean(ni1.data(:,:,:,5:end),4); 
        % this worked best for subject 2, scan 2

        ni1.dim=ni1.dim(1:3);
        ni1.pixdim=ni1.pixdim(1:3);
        niAnatomy.pixdim=niAnatomy.pixdim(1:3); % only uses the first three dimensions

        % allign functionals to the T1:
        acpcXform = dtiRawAlignToT1(ni1,niAnatomy,[], [], false, 1); % last 1 adds a figure
        % this saves the reallignment matrix in the folder of the functionals, 
    end
    
%     %%%%% coregister the venogram to the T1:
%     niAnatomy.pixdim=niAnatomy.pixdim(1:3); % only uses the first three
%     niVeno.pixdim=niVeno.pixdim(1:3); % only uses the first three
%     % allign Veno to the T1:
%     acpcXformVeno = dtiRawAlignToT1(niVeno,niAnatomy,[], [], [], 1); % last 1 adds a figure
%     % this saves the reallignment matrix in the folder of the venogram, 

end

% Only if the current acpcXform is good enough, safe for use
% Do this for the FA = 25, then coregister other functionals to this one
% ... only use the first visualization step from Kendrick's code to check...
acpcXform_new = acpcXform;
save(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']),'acpcXform_new')


%% Check the coregistration
s = 2;
s_info = bb_subs(s);
subj=s_info.subj;

% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));

scan_nr = 1; % selects the functional scan, note that we're using the second scan here, with FA = 25|20, the other FA do not coregister well to the T1
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};

fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
ni = niftiRead(fmri);
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']),'acpcXform_new')

curPos = [1,10,-30]; 
sliceThisDim = 2; 
% imDims=[-90 -120 -60; 90 130 90];
imDims=[-90 -120 -120; 90 130 90];

%%%% Overlay 1: functionals and anatomy
niFunc = ni;
niFunc.data = mean(ni.data(:,:,:,5:end),4); % overlay the mean functional
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform_new,sliceThisDim,imDims,curPos)





%% The rest seems to be unnecessary

%% NOTE, Kendrick's FUNCTIONS DO NOT WORK ANYMORE!
% there is somethin with allignvolumedata that makes my Matlab crash...
% specifically during alignvolumedata_auto

%%
%% refine acpcXform by using Kendrick's alignment code:
%%
%% get data

clear all
close all
% dDir = '/biac4/wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

s_nr = 2;
scan_nr = 2;
  
s_info = bb_subs(s_nr);
subj=s_info.subj;

% Get the functionals
scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};
ni = niftiRead(fullfile(dDir,subj,scan, [scanName '.nii.gz']));

% use the first functional
ni.data = ni.data(:,:,:,1);

% Get the anatomical:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii']));
% % use other functional instead of anatomical:
% niAnatomy = niftiRead(fullfile(dDir,subj,s_info.scan{2}, [s_info.scanName{2} '.nii.gz']));
% niAnatomy.data = niAnatomy.data(:,:,:,1);

% get the initial coregistration matrix
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new;

%% Open alignvolumedata gui from Kendrick
% this still works

% modify initial alignment to KNK format
% rxAlignment = niAnatomy.qto_ijk * ni.qto_xyz; % for FA48 ?
rxAlignment = niAnatomy.qto_ijk * acpcXform; % start out with best SPM alignment

knk.TORIG = rxAlignment;
knk.trORIG = matrixtotransformation(knk.TORIG,0,niAnatomy.pixdim(1:3),size(ni.data),size(ni.data) .* ni.pixdim(1:3));

volpre = ni.data(:,:,:,1);
refpre = niAnatomy.data;

% Open alignment gui
alignvolumedata(refpre,niAnatomy.pixdim(1:3),volpre,ni.pixdim(1:3),knk.trORIG);

%%
%% Allign functional to a good functional with SPM
%%
clear all
close all
% dDir = '/biac4/wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

s_nr = 2;
scan_nr = 2;
ref_scan_nr = 2;

s_info = bb_subs(s_nr);
subj=s_info.subj;

% Get the new ref scan:
%     niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.scan{ref_scan_nr}, [s_info.scanName{ref_scan_nr} '.nii.gz']));
niAnatomy.data = niAnatomy.data(:,:,:,1);

%%%%% coregister the functionals to the ref funx:

scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};

fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
ni = niftiRead(fmri);

%%%%% use the first nifti to allign, this one has the most structural info:
ni1=ni;
ni1.data=ni1.data(:,:,:,1);
ni1.dim=ni1.dim(1:3);
ni1.pixdim=ni1.pixdim(1:3);
niAnatomy.pixdim=niAnatomy.pixdim(1:3); % only uses the first three dimensions

% allign functionals to the T1:
acpcXform = dtiRawAlignToT1(ni1,niAnatomy,[], [], false, 1); % last 1 adds a figure
% this saves the reallignment matrix in the folder of the functionals, 

% load ref scan acpc x-form
ref_acpc = load(fullfile(dDir,subj,s_info.scan{ref_scan_nr},...
    [s_info.scanName{ref_scan_nr} 'AcpcXform_new.mat']));

% now fix the acpcXform such that it goes to T1 space
acpcXform_new = ref_acpc.acpcXform_new * niAnatomy.qto_ijk * acpcXform; % funx -> ref xyz -> ref ijk -> ref acpc
save(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']),'acpcXform_new')

%%
%% NOW MAKE SOME FIGURES
%%

%% The T2* data are here.  
% dDir = '/biac4/wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

in_data = 'PPG';
s = 2;
scan_nr = 1;

curPos = [10,1,1]; 
sliceThisDim = 1; 
imDims=[-90 -120 -60; 90 130 90];
% imDims=[-90 -120 -120; 90 130 90];

s_info = bb_subs(s);
subj = s_info.subj;
    
% Get the anatomicals:
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii']));

% % Get the MRVenogram:
% niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));

scan=s_info.scan{scan_nr};
scanName=s_info.scanName{scan_nr};

fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
ni = niftiRead(fmri);

% load coregistration matrix (for the functionals):
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% load time series and associated time
ppgTSname = fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponse.nii.gz']);
ppgTS = niftiRead(ppgTSname); % ppg triggered time series
load(fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponseT.mat']),'t');

% load correlation between even and odd scans for colors of timeseries
ppgRname = fullfile(dDir,subj,scan,[scanName '_cod' in_data '.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

%%%% Overlay 1: functionals and anatomy
niFunc = ni;
niFunc.data = mean(ni.data(:,:,:,5:end),4); % overlay the mean functional
bbOverlayFuncAnat(niFunc,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)

%%%% Overlay 2: timeseries and anatomy
ppgTSplot = ppgTS;
if isequal(in_data,'PPG')
    ppgTSplot.data(:,:,:,t<-.1 | t>1)=[]; % plot these times from curve
elseif isequal(in_data,'RESP')
    ppgTSplot.data(:,:,:,t<-.2 | t>4)=[]; % plot these times from curve
end
% Scale time series amplitude by R, plots are generated with respect to the maximum.
niColor = ppgR; % use R2 map to scale later: r2_scale=sqrt(imgSlice2.^2);
maxTS = max(abs(ppgTSplot.data),[],4); % get the max of each curve
ppgTSplot.data = bsxfun(@rdivide,ppgTSplot.data,maxTS); % devide by the max of each curve (sets all curves to 1 max)
ppgTSplot.data = bsxfun(@times,ppgTSplot.data,niColor.data); % multiply by r^2 to set less reliable curves to zero

bbOverlayTimeseriesAnat(ppgTSplot,niColor,niAnatomy,acpcXform,sliceThisDim,imDims,curPos)

clear niColor ppgTSplot

%% plot MRV and timeseries
% Get the MRVenogram:
niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));
% load coregistration matrix (for the venogram):
xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));

bbOverlayTimeseriesVeno(ppgTSplot,niColor,niVeno,acpcXform,xf_veno.acpcXform,sliceThisDim,imDims,curPos)


%% plot MRV and T1

% Get the MRVenogram:
niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));
% load coregistration matrix (for the venogram):
xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));

% then use dtiGetSlice to get the same slice from 2 sets
curPos = [-1,10,1]; 
sliceThisDim = 1; 
% imDims=[-90 -120 -60; 90 130 90];
imDims=[-90 -120 -120; 90 130 90];

% get a slice from the MRV
imgVol = niVeno.data;%ni.data(:,:,:,1);
img2std = xf_veno.acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [1 1 1];
[imgSlice1,x1,y1,z1]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
x1=x1(1,:)';
y1=y1(:,1);
z1=z1(1,:)';

% get the same slice for the anatomy
imgVol = niAnatomy.data;
img2std = niAnatomy.qto_xyz;
sliceNum =curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [1 1 1];
[imgSlice,x,y,z]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
x=x(1,:)';
y=y(:,1);
z=z(1,:)';

% for x and y for plotting:
if sliceThisDim==1
    x1=z1; x=z;
    % flip x and y
    x1_t=x1;
    x1=y1; 
    y1=x1_t;
    x_t=x;
    x=y; 
    y=x_t;
    
    % and for the images
    imgSlice1=imrotate(imgSlice1,90);
    imgSlice=imrotate(imgSlice,90);
elseif sliceThisDim==2
    y1=x1; y=x;
    x1=z1; x=z;
elseif sliceThisDim==3
    % x and y stay the same
end

% now plot stuff:
f=figure;
cm=colormap(jet);
close(f)

figure('Position',[0 0 1200 600])
subplot(1,2,1)
imagesc(x1,y1,imgSlice1,[0 1500])
hold on
axis image

subplot(1,2,2)
% show the background:
image(x,y,cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:)));
hold on
% colormap gray
axis image

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,subj,scan,[scanName '_Veno_T1_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))
% print('-painters','-r300','-depsc',fullfile(dDir,subj,scan,[scanName '_Veno_T1_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))


%% plot first and second functional and T1

% then use dtiGetSlice to get the same slice from 2 sets
% curPos = [1,10,-20]; 
curPos = [1,10,-20]; 
sliceThisDim = 1; 
% imDims=[-90 -120 -60; 90 130 90];
imDims=[-90 -120 -120; 90 130 90];

% get a slice from the functionals
imgVol1 = ni.data(:,:,:,1);
imgVol1 = imgVol1/max(imgVol1(:)); % normalize maximum to 1
imgVol2 = ni.data(:,:,:,2);
imgVol2 = imgVol2/max(imgVol2(:)); % normalize maximum to 1
img2std = acpcXform;
sliceNum = curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [4 4 4];
[imgSlice1,x1,y1,z1]=dtiGetSlice(img2std, imgVol1, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
[imgSlice2]=dtiGetSlice(img2std, imgVol2, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
x1=x1(1,:)';
y1=y1(:,1);
z1=z1(1,:)';

% get the same slice for the anatomy
imgVol = niAnatomy.data;
img2std = niAnatomy.qto_xyz;
sliceNum =curPos(sliceThisDim);
interpType = 'n';
mmPerVox = [1 1 1];
[imgSlice,x,y,z]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
x=x(1,:)';
y=y(:,1);
z=z(1,:)';

% for x and y for plotting:
if sliceThisDim==1
    x1=z1; x=z;
    % flip x and y
    x1_t=x1;
    x1=y1; 
    y1=x1_t;
    x_t=x;
    x=y; 
    y=x_t;
    
    % and for the images
    imgSlice1=imrotate(imgSlice1,90);
    imgSlice2=imrotate(imgSlice2,90);
    imgSlice=imrotate(imgSlice,90);
elseif sliceThisDim==2
    y1=x1; y=x;
    x1=z1; x=z;
elseif sliceThisDim==3
    % x and y stay the same
end

% now plot stuff:
f=figure;
cm=colormap(jet);
close(f)

figure('Position',[0 0 1200 600])
subplot(1,3,1)
imagesc(x1,y1,imgSlice1)
hold on
axis image
title('volume 1')

subplot(1,3,2)
imagesc(x1,y1,imgSlice2)
hold on
axis image
title('volume 2')

subplot(1,3,3)
image(x,y,cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:)));
hold on
colormap gray
axis image
title('T1')

set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',fullfile(dDir,subj,scan,[scanName '_FunxandT1_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))





