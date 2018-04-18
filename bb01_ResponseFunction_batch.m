clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses


%% Base data directory 

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';


%% Preprocess the data: PPG time series and correlation

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

for s = 4
    s_info = bb_subs(s);
    subj = s_info.subj;
    for scan_nr = [2 4 5 7 8]%:length(s_info.scan)
        scan = s_info.scan{scan_nr};
        scanName = s_info.scanName{scan_nr};

        fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
        ni = niftiRead(fmri);

        % compute the PPG triggered response matrix for the whole brain and save as a nifti:
        [response_matrix,t,response_matrix_odd,response_matrix_even,response_matrix_sterr] ...
            = bbResponse2physio(ni);

        % safe time T:
        save(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT']),'t')

        % save average of all heartbeats:
        ni1 = ni;
        ni1.data = response_matrix;
        ni1.fname = [scanName '_PPGtrigResponse.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1

        % save standard error of all heartbeats:
        ni1 = ni;
        ni1.data = response_matrix_sterr;
        ni1.fname = [scanName '_PPGtrigResponse_std.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
        
        % save average of all odd heartbeats:
        ni1 = ni;
        ni1.data = response_matrix_odd;
        ni1.fname = [scanName '_PPGtrigResponse_odd.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
        
        % save average of all even heartbeats:
        ni1 = ni;
        ni1.data = response_matrix_even;
        ni1.fname = [scanName '_PPGtrigResponse_even.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
    end
    
    for scan_nr = [1:9]%1:length(s_info.scan)
        scan = s_info.scan{scan_nr};
        scanName = s_info.scanName{scan_nr};

        ni_odd = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_odd.nii.gz']));
        ni_even = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_even.nii.gz']));
        
        % compute the correlation between the even and odd repetitions of
        % the PPG triggered BOLD signal for the whole brain:
        [out_r_map] = bbCorrelate2physio(ni_odd,ni_even);
        % save as nifti
        ni1 = ni_odd;
        ni1.data = out_r_map;
        ni1.fname = [scanName '_corrPPG.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1 out_r_map

        % compute the cod between the even and odd repetitions of
        % the PPG triggered BOLD signal for the whole brain:
        [out_R_map] = bbCod2physio(ni_odd,ni_even);
        % save as nifti
        ni2 = ni_odd;
        ni2.data = out_R_map;
        ni2.fname = [scanName '_codPPG.nii.gz'];
        niftiWrite(ni2,fullfile(dDir,subj,scan,ni2.fname))
        clear ni2 out_R_map
       
    end
end

%%
%% Preprocess the data: RESP time series and correlation

% The pixdim field in the ni structure has four dimensions, three spatial
% and the fourth is time in seconds.

for s=2
    s_info = bb_subs(s);
    subj=s_info.subj;
    for scan_nr=2%1:length(s_info.scan)
        scan=s_info.scan{scan_nr};
        scanName=s_info.scanName{scan_nr};

        fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
        ni = niftiRead(fmri);

        % compute the PPG triggered response matrix for the whole brain and save as a nifti:
        [response_matrix,t,response_matrix_odd,response_matrix_even] = bbResponse2physio(ni,[],[-0.5 5],'resp');

        % safe time T:
        save(fullfile(dDir,subj,scan,[scanName '_RESPtrigResponseT']),'t')

        % save average of all heartbeats:
        ni1=ni;
        ni1.data=response_matrix;
        ni1.fname=[scanName '_RESPtrigResponse'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
        
        % save average of all odd heartbeats:
        ni1=ni;
        ni1.data=response_matrix_odd;
        ni1.fname=[scanName '_RESPtrigResponse_odd'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
        
        % save average of all even heartbeats:
        ni1=ni;
        ni1.data=response_matrix_even;
        ni1.fname=[scanName '_RESPtrigResponse_even'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1

    end
    for scan_nr=2%1:length(s_info.scan)
        scan=s_info.scan{scan_nr};
        scanName=s_info.scanName{scan_nr};

        ni_odd = niftiRead(fullfile(dDir,subj,scan,[scanName '_RESPtrigResponse_odd']));
        ni_even = niftiRead(fullfile(dDir,subj,scan,[scanName '_RESPtrigResponse_even']));
        
        % compute the correlation/ reliability with the RESP for the whole brain and save as a nifti:
        [out_r_map]=bbCorrelate2physio(ni_odd,ni_even);

        % save as nifti
        ni1=ni_odd;
        ni1.data=out_r_map;
        ni1.fname=[scanName '_corrRESP.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
    end
end

%% Whole brain measurement of reliability

% Select a subject and scan nummer
s_nr = 4;
scan_nr = 3;

roisToSegment = {[...
    4 % left lateral ventricle
    5 % left inferior lateral ventricle
    14 % 3rd ventricle
    15 % 4th ventricle
    24 % CSF
    31 % left choroid plexus
    43 % right lateral ventricle
    44 % right inferior lateral ventricle
    63 % right choroid plexus
    72],... % 5th ventricle
    [2 % left white matter
    41],... % right white matter
    [3 % left gray matter
    42]}; % right gray matter

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

% load coregistration matrix for the functionals:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% Get the MRVenogram:
niVeno = niftiRead(fullfile(dDir,subj,subs.veno,[subs.venoName '.nii.gz']));
% load coregistration matrix (for the venogram):
xf_veno=load(fullfile(dDir,subj,subs.veno,[subs.venoName 'AcpcXform.mat']));

% Anatomical
anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii.gz']);
niAnatomy = niftiRead(anat);

% COD even/odd
ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

% Segmentation file:
segName = fullfile(dDir,subj,scan,[scanName '_r_aseg_auto.nii.gz']);
niSeg = niftiRead(segName);
 
% Veno file:
venoName = fullfile(dDir,subj,scan,[scanName '_r_veno.nii.gz']);
niVeno = niftiRead(venoName);

% % Get voxels with mean signal above brain_th
% brain_th = 50;
% brain_vect = mean(ni.data(:,:,:,4:end),4);
% brain_vect = brain_vect>brain_th;

clear out

totalNrVox = 0;
for rr = 1:length(roisToSegment)
    niROIcode = roisToSegment{rr};
    voxels_inroi = ismember(niSeg.data(:),niROIcode);
    out(rr).roi_ppgR = ppgR.data(voxels_inroi);
    totalNrVox = totalNrVox+length(find(voxels_inroi>0));
end

voxels_inroi = niVeno.data(:)>500;
out(rr+1).roi_ppgR = ppgR.data(voxels_inroi);

r_th = .3;

figure('Position',[0 0 500 400])
hold on
for rr = 1:length(out)
    bar(rr,length(find(out(rr).roi_ppgR>r_th))./length(out(rr).roi_ppgR));
end





%% OLD:
%% Whole brain measurement of reliability

%%% TODO: change this to include code from bb02_CreateROIs to get a better
%%% estimate for voxels in each roi

% Select a subject and scan nummer
s_nr = 4;
scan_nr = 9;

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

% load coregistration matrix for the functionals:
load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
acpcXform = acpcXform_new; clear acpcXform_new

% Get the MRVenogram:
niVeno = niftiRead(fullfile(dDir,subj,subs.veno,[subs.venoName '.nii.gz']));
% load coregistration matrix (for the venogram):
xf_veno=load(fullfile(dDir,subj,subs.veno,[subs.venoName 'AcpcXform.mat']));

% Anatomical
anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii.gz']);
niAnatomy = niftiRead(anat);

% COD even/odd
ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
ppgR = niftiRead(ppgRname); % correlation with PPG

% Get voxels with mean signal above brain_th
brain_th = 50;
brain_vect = mean(ni.data(:,:,:,4:end),4);
brain_vect = brain_vect>brain_th;

% Which percentage of voxels with mean signal above threshold had
% replicable brain-beat curves?
total_nr = length(find(brain_vect==1));

r_th = .3;
reliable_nr = length(find(brain_vect == 1 & ppgR.data>r_th));
roisToSegmentNames = {...
    'L_lat_ventr'
    'L_inferior_lat_ventr'
    '3rd_ventr'
    '4th_ventr'
    'CSF'
    'L_choroid_plexus'
    'R_lat_ventr'
    'R_inferior_lat_ventr'
    'R_choroid_plexus'
    '5th_ventricle'
    'L_white'
    'L_gray'
    'R_white'
    'R_gray'};

clear out

for rr = 1:length(roisToSegmentNames)

    niROIname = roisToSegmentNames{rr};
    niROI = niftiRead(fullfile(dDir,subj,'freesurfer','nii',[niROIname '.nii.gz']));

    % get ROI indices:
    [xx,yy,zz] = ind2sub(size(niROI.data),find(niROI.data>0));

    % now ROI indices to ACPC (mm):
    xyz_acpc = mrAnatXformCoords(niROI.qto_xyz, [xx,yy,zz]);
    clear xx yy zz % housekeeping

    % now ACPC coordinates to functional indices:
    ijk_func = mrAnatXformCoords(inv(acpcXform), xyz_acpc);
    ijk_func = round(ijk_func); % round to get indices
    ijk_func = unique(ijk_func,'rows'); % only take unique voxels

    %%%% check for coordinates in functional space
    z_slice = [5 10 15 18 20 23 25 26 27 30 33 36];
    figure
    for kk = 1:length(z_slice)       
        subplot(3,4,kk)
        imagesc(ni.data(:,:,z_slice(kk))),hold on
        xyz_plot = ijk_func(ijk_func(:,3)==z_slice(kk),:);
        plot(xyz_plot(:,2),xyz_plot(:,1),'r.')
        axis image 
    end

    % remove ijk smaller than 1 - indices do not fit in fMRI data
    ijk_func(sum(ijk_func<1,2)==1,:) = [];

    out(rr).roi_ppgR = zeros(size(ijk_func,1),1);
    for kk = 1:size(ijk_func,1)
        out(rr).roi_ppgR(kk) = ppgR.data(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3));
    end
end

% Get venogram ROI
niVeno.data = single(niVeno.data);
% get ROI indices:
[xx,yy,zz] = ind2sub(size(niVeno.data),find(niVeno.data>1500));

% now ROI indices to ACPC (mm):
xyz_acpc = mrAnatXformCoords(niVeno.qto_xyz, [xx,yy,zz]);
clear xx yy zz % housekeeping

% now ACPC coordinates to functional indices:
ijk_func = mrAnatXformCoords(inv(acpcXform), xyz_acpc);
ijk_func = round(ijk_func); % round to get indices
ijk_func = unique(ijk_func,'rows'); % only take unique voxels
%%% check for coordinates in functional space
% z_slice = [5 10 15 18 20 23 25 26 27 30 33 36];
% figure
% for kk = 1:length(z_slice)       
%     subplot(3,4,kk)
%     imagesc(ni.data(:,:,z_slice(kk))),hold on
%     xyz_plot = ijk_func(ijk_func(:,3)==z_slice(kk),:);
%     plot(xyz_plot(:,2),xyz_plot(:,1),'r.')
%     axis image 
% end

% remove ijk smaller than 1 - indices do not fit in fMRI data
ijk_func(sum(ijk_func<1,2)==1,:) = [];

out(rr+1).roi_ppgR = zeros(size(ijk_func,1),1);
for kk = 1:size(ijk_func,1)
    out(rr+1).roi_ppgR(kk) = ppgR.data(ijk_func(kk,1),ijk_func(kk,2),ijk_func(kk,3));
end

ROInames = {...
    'lLV'
    'lILV'
    '3V'
    '4V'
    'CSF'
    'lCh'
    'rLV'
    'rILV'
    'rCh'
    '5V'
    'lW'
    'lG'
    'rW'
    'rG'
    'veno'};

figure('Position',[0 0 500 400])
hold on
for rr = 1:length(ROInames)
    bar(rr,length(find(out(rr).roi_ppgR>r_th))./length(out(rr).roi_ppgR));
end

% Add total
bar(rr+1,reliable_nr./total_nr)

set(gca,'XTick',[1:rr+1],'XTickLabel',{ROInames{:},'all'})
ylim([0 1])
ylabel(['% voxels with R>' num2str(r_th,3)])
title(['flip' int2str(subs.scanFA{scan_nr})])

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',[dDir './figures/reliable/subj' int2str(s_nr) '_scan' int2str(scan_nr)])
print('-painters','-r300','-depsc',[dDir './figures/reliable/subj' int2str(s_nr) '_scan' int2str(scan_nr)])
