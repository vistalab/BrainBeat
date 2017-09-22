clear all
close all

% script produced transformation matrix between the anatomy and functional
% scans, walk through step-by-step

addpath(genpath('~/Documents/m-files/kendrick/'))

%% Base data directory on a Mac mounting biac4 (wandell's machine)
% dDir = '/Volumes/biac4-wandell/data/BrainBeat/data';
% dDir = '/biac4/wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% chdir(dDir)

%% coregistration functionals to T1 (course)
% this can work well for FA - 25, but not for FA = 36/48, thus: allign FA =
% 25, and align 36/48 to this one with code below
 
for s = 3
    s_info = bb_subs(s);
    subj=s_info.subj;

    % Get the anatomicals:
    niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii']));

%     % Get the MRVenogram:
%     niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii']));

    %%%%% coregister the functionals to the T1:

    for scan_nr = 3
        scan=s_info.scan{scan_nr};
        scanName=s_info.scanName{scan_nr};

        fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
        ni = niftiRead(fmri);

        ni1=ni;
        %%%%% use the first nifti to allign, this one has the most structural info:
        ni1.data=ni1.data(:,:,:,1);
        %%%%% or try the mean of all nifti's for the FA48:
%         ni1.data=mean(ni1.data(:,:,:,4:end),4);

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
% ... only use sthe first visualization step from Kendick's code to check...
% acpcXform_new = acpcXform;
% save(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']),'acpcXform_new')


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
scan_nr = 3;
  
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

%% Open alignvolumedata gui
% this still works

% modify initial alignment to KNK format
% rxAlignment = niAnatomy.qto_ijk * ni.qto_xyz; % for FA48 ?
rxAlignment = niAnatomy.qto_ijk * acpcXform; % start out with best SPM alignment

knk.TORIG = rxAlignment;
knk.trORIG = matrixtotransformation(knk.TORIG,0,niAnatomy.pixdim(1:3),size(ni.data),size(ni.data) .* ni.pixdim(1:3));

% test whether this helps, did not in one case (?):
% volpre = preconditionvolume(ni.data(:,:,:,1),[],[],[99 1/3]);
% volpre = preconditionvolume(ni.data(:,:,:,1));
% refpre = preconditionvolume(niAnatomy.data(:,:,:));
% close all
volpre = ni.data(:,:,:,1);
refpre = niAnatomy.data;

% Open alignment gui
alignvolumedata(refpre,niAnatomy.pixdim(1:3),volpre,ni.pixdim(1:3),knk.trORIG);

%% Define ellipse
[~,mn,sd] = defineellipse3d(double(volpre));
% s_nr = 2, scan_nr = 1 
% mn = [0.5048    0.5059    0.7692];
% sd = [0.2878    0.3105    0.0714];
mn = [0.5048    0.5589    0.6512];
sd = [0.2524    0.2782    0.2189];

%% Automatic alignments (coarse)
% DEFECT
useMI = false;
alignvolumedata_auto(mn,sd,0,[4 4 2],[],[],[],useMI) % rigid body, course, mutual information metric

%% Automatic alignments (fine)
% DEFECT
alignvolumedata_auto(mn,sd,0,[1 1 1],[],[],[],useMI) % rigid body, course, mutual information metric

tr = alignvolumedata_exporttransformation;
% to check again:
% alignvolumedata(refpre,niAnatomy.pixdim(1:3),volpre,ni.pixdim(1:3),tr);

%%

T = transformationtomatrix(tr,0,niAnatomy.pixdim(1:3));

% get T back to funx only:
% new acpcXform:
acpcXform_new = niAnatomy.qto_xyz * T;

% save new coregistration matrix
save(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']),'acpcXform_new')

%%
%% END TO NOTE, Kendrick's FUNCTIONS DO NOT WORK ANYMORE!
%%


%%
%% Allign functional to a good functional with SPM
%%
clear all
close all
% dDir = '/biac4/wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

s_nr = 3;
scan_nr = 5;
ref_scan_nr = 2;

s_info = bb_subs(s_nr);
subj=s_info.subj;

% Get the new ref scan:
%     niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii.gz']));
niAnatomy = niftiRead(fullfile(dDir,subj,s_info.scan{ref_scan_nr}, [s_info.scanName{ref_scan_nr} '.nii']));
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
ref_acpc = load(fullfile(dDir,subj,s_info.scan{ref_scan_nr},[s_info.scanName{ref_scan_nr} 'AcpcXform_new.mat']));

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

make_fig = [0 1]; % 0 do not, 1 do make, 2 also save figures: 
% [0 2] makes and saves figure 2 

in_data = 'PPG';

for s = 3%:3;

    s_info = bb_subs(s);
    subj = s_info.subj;
    
    % Get the anatomicals:
    niAnatomy = niftiRead(fullfile(dDir,subj,s_info.anat,[s_info.anatName '.nii']));

    % % Get the MRVenogram:
    % niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii.gz']));

for scan_nr = 3% 1:length(s_info.scan)%;
    scan=s_info.scan{scan_nr};
    scanName=s_info.scanName{scan_nr};

    fmri = fullfile(dDir,subj,scan, [scanName '.nii.gz']);
    ni = niftiRead(fmri);

    % Load the correlation with heartbeat (made with bbCorrelate2physio):
    ppgRname = fullfile(dDir,subj,scan,[scanName '_corr' in_data '.nii.gz']);
    ppgR = niftiRead(ppgRname); % correlation with PPG

    % Load the timeseries around PPG peak (made with bbResponse2physio):
    ppgTSname = fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponse.nii']);
    ppgTS = niftiRead(ppgTSname); % ppg triggered time series
    load(fullfile(dDir,subj,scan,[scanName '_' in_data 'trigResponseT.mat']),'t');

    % load coregistration matrix (for the functionals):
    load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
    acpcXform = acpcXform_new;
    
    % % load coregistration matrix (for the venogram):
    % xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));

    %%%% now overlay r-map with anatomy

    % then use dtiGetSlice to get the same slice from 2 sets
    sliceThisDim = 1; 
    
    if s==2
        imDims = [-90 -120 -120; 90 130 90]; 
        curPos = [1,10,-20]; 
        curPos = [-29,-14,-50]; % 03
        curPos = [-3,30,-43]; % 03
    elseif s==3
        imDims = [-90 -120 -100; 90 130 110]; 
        curPos = [1,4,38]; 
    end
    
    % functionals to ACPC space
    % settings:
    img2std = acpcXform;
    sliceNum = curPos(sliceThisDim);
    interpType = 'n';
    mmPerVox = [4 4 4];

    % functionals to ACPC
    imgVol = ni.data(:,:,:,1); % do the first
    imgVol = imgVol/max(imgVol(:));
    [imgSlice1,x1,y1,z1]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
    if sliceThisDim == 1 || sliceThisDim == 3
        x1=x1(1,:)';
        y1=y1(:,1);
    elseif sliceThisDim == 2
        x1=x1(:,1);
        y1=y1(1,:)';
    end
    z1=z1(1,:)';

    % correlation map to ACPC
    imgVol = ppgR.data;
    [imgSlice2]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);

    % PPG timeSeries to ACPC
    imgVol = ppgTS.data;
    [imgSlice3]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);

    % Anatomy to ACPC
    imgVol = niAnatomy.data;
    img2std = niAnatomy.qto_xyz;
    sliceNum =curPos(sliceThisDim);
    interpType = 'n';
    mmPerVox = [1 1 1];
    [imgSlice,x,y,z]=dtiGetSlice(img2std, imgVol, sliceThisDim, sliceNum,imDims,interpType, mmPerVox);
    if sliceThisDim == 1 || sliceThisDim == 3
        x=x(1,:)';
        y=y(:,1);
    elseif sliceThisDim == 2
        x=x(:,1);
        y=y(1,:)';
    end
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
        imgSlice3=imrotate(imgSlice3,90);
        imgSlice=imrotate(imgSlice,90);
    elseif sliceThisDim==2
        y1=z1; y=z;

        % rotate the images
        imgSlice1=imrotate(imgSlice1,90);
        imgSlice2=imrotate(imgSlice2,90);
        imgSlice3=imrotate(imgSlice3,90);
        imgSlice=imrotate(imgSlice,90);

    elseif sliceThisDim==3
        % x and y stay the same
    end

    % now plot stuff:
    f=figure;
    cm=colormap(jet);
    close(f)

    if make_fig(1) > 0
        figure('Position',[0 0 1200 1200])
        subplot(2,2,1)
        imagesc(x1,y1,imgSlice1)
        hold on
        axis image
        title('functional 1')
        
        subplot(2,2,2)
        image(x,y,cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:))); % background
        hold on
        axis image
        % tranform imgSlice1 (functionals) to colormap values that I want to use
        imgSlice1_color=cat(3,zeros(size(imgSlice1)),zeros(size(imgSlice1)),zeros(size(imgSlice1)));
        for k_x = 1:length(x1)
            for k_y = 1:length(y1)
                imgSlice1_color(k_y,k_x,:) = cm(1+floor(abs(imgSlice1(k_y,k_x))*64),:);
            end
        end
        h=image(x1,y1,imgSlice1_color); % overlay colormap on image
        set(h,'AlphaData',.2*ones(size(imgSlice1))) % make transparent
        
        subplot(2,2,3)
        imagesc(x1,y1,imgSlice2,[0 1])
        hold on
        axis image
        title('correlation map')
        
        subplot(2,2,4)
        image(x,y,cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:))); % background
        hold on
        axis image
        % tranform imgSlice2 (r-map) to colormap values that I want to use
        imgSlice1_color=cat(3,zeros(size(imgSlice2)),zeros(size(imgSlice2)),zeros(size(imgSlice2)));
        for k_x = 1:length(x1)
            for k_y = 1:length(y1)
                imgSlice1_color(k_y,k_x,:) = cm(1+floor(abs(imgSlice2(k_y,k_x))*64),:);
            end
        end
        h=image(x1,y1,imgSlice1_color); % overlay colormap on image
        set(h,'AlphaData',.2*ones(size(imgSlice1))) % make transparent
        
        if make_fig(1) == 2
            % alpha values do not save in png...
            set(gcf,'PaperPositionMode','auto')
        %     print('-painters','-r300','-dpng',fullfile(dDir,subj,scan,[scanName '_corrPPG_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))
%             print('-painters','-r300','-dpng',['./figures/test_' in_data '/' in_data 'corr_' subj '_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])
        end
    end
    
    if make_fig(2) > 0
    %%%%% overlay with time series
    f=figure;
    cm=colormap(jet);
    close(f)

    figure('Position',[0 0 1200 600])

    subplot(1,2,1)
    % show the background:
    imagesc(x1,y1,imgSlice2,[0 1]);
    hold on
    colormap gray
    axis image

    subplot(1,2,2)
    % show the background:
%     image(x,y,cat(3,imgSlice,imgSlice,imgSlice)/max(imgSlice(:)));
    imagesc(x,y,imgSlice/max(imgSlice(:)));
    set(gca,'CLim',[0 .3])
    hold on
    axis image
    
    if isequal(in_data,'PPG')
        a=squeeze(imgSlice3(:,:,t>-.2 & t<1)); % 3rd dimension is time
    elseif isequal(in_data,'RESP')
        a=squeeze(imgSlice3(:,:,t>-.2 & t<4)); % 3rd dimension is time
    end
    r2_scale=sqrt(imgSlice2.^2);

    for sub_p=1:2
        subplot(1,2,sub_p)
    for k=1:size(a,1)
        for m=1:size(a,2)
            r_curve=squeeze(a(k,m,:));
            
            % get plotting location
            k_x = y1(k);
            m_y = x1(m);

            % scale by strength of correlation with heartbeat:
            r_curve = r_curve/max(abs(r_curve(:))); % set max to 1
            r_curve = r_curve*r2_scale(k,m); % scale with r
            
            % scaling factor to make the curve fit in a 4 mm voxels 
            s_f=3; 

            c_use=cm(1+floor(r2_scale(k,m)*64),:);
            plot(m_y+s_f*[0:size(a,3)-1]/size(a,3)-2,k_x+s_f*r_curve,'Color',c_use)

    %         plot(m_y-2,k_x,'.','Color',[.5 .5 .5])
        end
    end
    end

    if make_fig(2) == 2
        set(gcf,'PaperPositionMode','auto')
    %     print('-painters','-r300','-dpng',fullfile(dDir,subj,scan,[scanName '_BBcurves_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))
        % print('-painters','-r300','-depsc',fullfile(dDir,subj,scan,[scanName '_BBcurves_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))
%         print('-painters','-r300','-dpng',['./figures/test_' in_data '/' in_data 'curves_' subj '_FA' int2str(s_info.scanFA{scan_nr}) '_scan' int2str(scan_nr) '_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))])
    end
    end
end
end

%% plot MRV and T1

% then use dtiGetSlice to get the same slice from 2 sets
curPos = [20,10,1]; 
sliceThisDim = 3; 
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

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,subj,scan,[scanName '_Veno_T1_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))
print('-painters','-r300','-depsc',fullfile(dDir,subj,scan,[scanName '_Veno_T1_view' int2str(sliceThisDim) '_slice' int2str(curPos(sliceThisDim))]))


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





