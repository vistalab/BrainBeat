

[~,dDir] = bbPath;

%
% ROIs were generated before BIDS conversion
%

%% Convert entire Freesurfer aseg.auto.mgz data to nifti in t1 space:

% clear all
% close all
% 

s = 6;
s_info = bb_subs(s);
subj = s_info.subj;

resample_type= 'weighted';

alignTo = fullfile(dDir,s_info.subj,s_info.anat,[s_info.anatName '.nii.gz']);
segmentFile = fullfile(dDir,s_info.subj,'freesurfer','mri','aseg.auto.mgz');
outfile = fullfile(dDir,s_info.subj,'freesurfer','nii','aseg.auto.nii.gz');

str = sprintf('!mri_convert  --out_orientation RAS --reslice_like %s -rt %s %s %s', alignTo, resample_type, segmentFile, outfile);

eval(str)

% read in the nifti
% ni = niftiRead(outfile);

%% Write aseg.auto in the space of a functional scan

% clear all
% close all

s = 6;

for scan_nr = [4 6]

    s_info = bb_subs(s);
    subj = s_info.subj;

    segm_file = fullfile(dDir,s_info.subj,'freesurfer','nii','aseg.auto.nii.gz');

    % read in the nifti
    ni = niftiRead(segm_file);

    % read the functional data
    scan = s_info.scan{scan_nr};
    scanName = s_info.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
    niFunc = niftiRead(fmri);

    % Load coregistration matrix:
    load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
    acpcXform = acpcXform_new; clear acpcXform_new

    % Create output file
    niFuncSeg = niFunc;
    niFuncSeg.data = niFuncSeg.data(:,:,:,1);
    niFuncSeg.data = zeros(size(niFuncSeg.data));
    niFuncSeg.fname = [scanName '_r_aseg_auto.nii.gz'];

    % for each func index, find a struct voxel
    for ii = 1:size(niFunc.data,1)
        for jj = 1:size(niFunc.data,2)
            for kk = 1:size(niFunc.data,3)
                % func voxel indices to acpc
                xyz_acpc = mrAnatXformCoords(acpcXform, [ii,jj,kk]);
                % acpc index to anat ijk
                ijk_anat = mrAnatXformCoords(ni.qto_ijk,xyz_acpc);
                % just take nearest neighbour value
                ijk_anat = round(ijk_anat);

                % check if the voxel index does not exceed anatomical
                if sum(ijk_anat<size(ni.data))==3 && sum(ijk_anat<1)==0
                    % get the image value and put it back in the functional
                    niFuncSeg.data(ii,jj,kk) = ni.data(ijk_anat(1),ijk_anat(2),ijk_anat(3));
                end
            end
        end
    end

    niftiWrite(niFuncSeg,fullfile(dDir,subj,scan,niFuncSeg.fname))
end

%% Write venogram in the space of a functional scan

% clear all
% close all


s = 6;

for scan_nr = [4 6]
    s_info = bb_subs(s);
    subj = s_info.subj;

    % Get the MRVenogram:
    niVeno = niftiRead(fullfile(dDir,subj,s_info.veno,[s_info.venoName '.nii.gz']));
    % load coregistration matrix (for the venogram):
    xf_veno=load(fullfile(dDir,subj,s_info.veno,[s_info.venoName 'AcpcXform.mat']));

    % read the functional data
    scan = s_info.scan{scan_nr};
    scanName = s_info.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
    niFunc = niftiRead(fmri);

    % Load coregistration matrix:
    load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
    acpcXform = acpcXform_new; clear acpcXform_new

    % Create output file
    niFuncVeno = niFunc;
    niFuncVeno.data = niFuncVeno.data(:,:,:,1);
    niFuncVeno.data = zeros(size(niFuncVeno.data));
    niFuncVeno.fname = [scanName '_r_veno.nii.gz'];

    % for each func index, find a struct voxel
    for ii = 1:size(niFunc.data,1)
        for jj = 1:size(niFunc.data,2)
            for kk = 1:size(niFunc.data,3)
                % func voxel indices to acpc
                xyz_acpc = mrAnatXformCoords(acpcXform, [ii,jj,kk]);
                % acpc index to anat ijk
                ijk_anat = mrAnatXformCoords(niVeno.qto_ijk,xyz_acpc);
                % just take neirest neighbour value
                ijk_anat = round(ijk_anat);

                % check if the voxel index does not exceed anatomical
                if sum(ijk_anat<size(niVeno.data))==3 && sum(ijk_anat<1)==0
                    % get the image value and put it back in the functional
                    niFuncVeno.data(ii,jj,kk) = niVeno.data(ijk_anat(1),ijk_anat(2),ijk_anat(3));
                end
            end
        end
    end

    niftiWrite(niFuncVeno,fullfile(dDir,subj,scan,niFuncVeno.fname))
end

%% Write SPM segmentation of the T1 in the space of a functional scan

% run SPM segment on anatomical

% clear all
% close all


s = 6;
s_info = bb_subs(s);
subj = s_info.subj;

for scan_nr = [4 6]
    % read gray, white and csf segmentations (respectively c1, c2, c3)
    segm_file1 = niftiRead(fullfile(dDir,s_info.subj,s_info.anat,['c1f' s_info.anatName '.nii']));
    segm_file2 = niftiRead(fullfile(dDir,s_info.subj,s_info.anat,['c2f' s_info.anatName '.nii']));
    segm_file3 = niftiRead(fullfile(dDir,s_info.subj,s_info.anat,['c3f' s_info.anatName '.nii']));
    segm_file1.data = double(segm_file1.data); segm_file2.data = double(segm_file2.data); segm_file3.data = double(segm_file3.data);     
    % set max back to 1 to get probabilities from 0 to 1
    segm_file1.data = segm_file1.data./max(segm_file1.data(:));
    segm_file2.data = segm_file2.data./max(segm_file2.data(:));
    segm_file3.data = segm_file3.data./max(segm_file3.data(:));
    
    % read the functional data
    scan = s_info.scan{scan_nr};
    scanName = s_info.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
    niFunc = niftiRead(fmri);

    % Load coregistration matrix:
    load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
    acpcXform = acpcXform_new; clear acpcXform_new

    % Create output file
    niFuncSeg = niFunc;
    niFuncSeg.data = niFuncSeg.data(:,:,:,1);
    niFuncSeg.data = zeros(size(niFuncSeg.data));
    niFuncSeg.fname = [scanName '_spmSeg.nii.gz'];
    
    % for each func index, find a struct voxel
    for ii = 1:size(niFunc.data,1)
        for jj = 1:size(niFunc.data,2)
            for kk = 1:size(niFunc.data,3)
                % func voxel indices to acpc
                xyz_acpc = mrAnatXformCoords(acpcXform, [ii,jj,kk]);
                % acpc index to anat ijk
                ijk_anat = mrAnatXformCoords(segm_file1.qto_ijk,xyz_acpc);
                % just take neirest neighbour value
                ijk_anat = round(ijk_anat);
                
                % check if the voxel index does not exceed anatomical
                if sum(ijk_anat<size(segm_file1.data))==3 && sum(ijk_anat<1)==0
                    % get the gray probability value
                    c1 = segm_file1.data(ijk_anat(1),ijk_anat(2),ijk_anat(3));
                    c2 = segm_file2.data(ijk_anat(1),ijk_anat(2),ijk_anat(3));
                    c3 = segm_file3.data(ijk_anat(1),ijk_anat(2),ijk_anat(3));
                    
                    if max([c1 c2 c3]) > .3 % only pick a tissue if decent probability
                        [~,max_pos] = max([c1 c2 c3]); % pick 1, 2, 3 for gray, white, csf
                        niFuncSeg.data(ii,jj,kk) = max_pos; % put it back in the functional space:
                    end
                end
            end
        end
    end

    niftiWrite(niFuncSeg,fullfile(dDir,subj,scan,niFuncSeg.fname))
end

%% Integrate all segmentations into 1 file for each scan: 

% clear all
% close all


s = 6;
s_info = bb_subs(s);
subj = s_info.subj;

for scan_nr = [4 6]%:3
    % read the functional data
    scan = s_info.scan{scan_nr};
    scanName = s_info.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
    niFunc = niftiRead(fmri);

    % freesurfer segmentation
    niFs = niftiRead(fullfile(dDir,subj,scan,[scanName '_r_aseg_auto.nii.gz']));
    
    % SPM segmentation
    niSPM = niftiRead(fullfile(dDir,subj,scan,[scanName '_spmSeg.nii.gz']));
    
    % venogram:
    niVeno = niftiRead(fullfile(dDir,subj,scan,[scanName '_r_veno.nii.gz']));
    
    % these are the freesurfer labels we need:
    FS_roisToSegment = {...
        [3 % left gray matter
        42],...% right gray matter 
        [2 % left white matter
        41],... % right white matter 
        [4 % left lateral ventricle
        5 % left inferior lateral ventricle
        14 % 3rd ventricle
        15 % 4th ventricle
        24 % CSF
        31 % left choroid plexus
        43 % right lateral ventricle
        44 % right inferior lateral ventricle
        63 % right choroid plexus
        72]... % 5th ventricle
        }; 
    
    % create output:
    newLabels = niSPM;
    newLabels.fname = [scanName '_combineSegm.nii.gz'];
    
    % remove GM and WM from SPM, get this from freesurfer:
    newLabels.data(newLabels.data==1 | newLabels.data==2) = 0;
    % give CSF label 4
    newLabels.data(newLabels.data==3) = 4;
    % get freesurfer labels and make gm = 1, wm = 2 and ventricles = 3
    for rr = 1:length(FS_roisToSegment)
        newLabels.data(ismember(niFs.data,FS_roisToSegment{rr})) = rr;
    end
    % give veins labels
    veno_th = 800;
    newLabels.data(niVeno.data>veno_th & niSPM.data>0) = 5;
    roiNames = {'GM','WM','Ventricles','CSF','Veno'};
    
    niftiWrite(newLabels,fullfile(dDir,subj,scan,newLabels.fname));
end

%%
%% Write DKT atlas in the space of a functional scan
%%
%% this is used for Figures 4 and 5

% clear all
% close all

% note that that aparc.DKTatlas+aseg.nii.gz is created through mri_convert

s = 6;

for scan_nr = [1]%[1:3]

    s_info = bb_subs(s);
    subj = s_info.subj;

    segm_file = fullfile(dDir,s_info.subj,'freesurfer','mri','aparc.DKTatlas+aseg.nii.gz');

    % read in the nifti
    ni = niftiRead(segm_file);

    % read the functional data
    scan = s_info.scan{scan_nr};
    scanName = s_info.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
    niFunc = niftiRead(fmri);

    % Load coregistration matrix:
    load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
    acpcXform = acpcXform_new; clear acpcXform_new

    % Create output file
    niFuncSeg = niFunc;
    niFuncSeg.data = niFuncSeg.data(:,:,:,1);
    niFuncSeg.data = zeros(size(niFuncSeg.data));
    niFuncSeg.fname = [scanName '_r_DKTatlas_aseg.nii.gz'];

    % for each func index, find a struct voxel
    for ii = 1:size(niFunc.data,1)
        for jj = 1:size(niFunc.data,2)
            for kk = 1:size(niFunc.data,3)
                % func voxel indices to acpc
                xyz_acpc = mrAnatXformCoords(acpcXform, [ii,jj,kk]);
                % acpc index to anat ijk
                ijk_anat = mrAnatXformCoords(ni.qto_ijk,xyz_acpc);
                % just take nearest neighbour value
                ijk_anat = round(ijk_anat);

                % check if the voxel index does not exceed anatomical
                if sum(ijk_anat<size(ni.data))==3 && sum(ijk_anat<1)==0
                    % get the image value and put it back in the functional
                    niFuncSeg.data(ii,jj,kk) = ni.data(ijk_anat(1),ijk_anat(2),ijk_anat(3));
                end
            end
        end
    end

    niftiWrite(niFuncSeg,fullfile(dDir,subj,scan,niFuncSeg.fname))
end
