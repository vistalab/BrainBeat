%% get ventricle ROI
clear all
close all

% dDir = '/biac4/wandell/data/BrainBeat/data';
dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

s = 3;
s_info = bb_subs(s);
subj = s_info.subj;

roisToSegment = [...
    4 % left lateral ventricle
    5 % left inferior lateral ventricle
    14 % 3rd ventricle
    15 % 4th ventricle
    24 % CSF
    31 % left choroid plexus
    43 % right lateral ventricle
    44 % right inferior lateral ventricle
    63 % right choroid plexus
    72]; % 5th ventricle

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
    '5th_ventricle'};

%%%% MAKE NII DIR IN FREESURFER FOLDER

%% make a nifti file from aseg.auto.mgz - then write only the ROIs all together in 1 file
resample_type = 'weighted';
alignTo = fullfile(dDir,s_info.subj,s_info.anat,[s_info.anatName '.nii']);
segmentFile = fullfile(dDir,s_info.subj,'freesurfer','mri','aseg.auto.mgz');
outfile = fullfile(dDir,s_info.subj,'freesurfer','nii','ventricles.nii.gz');

% make nifti from aseg
str = sprintf('!mri_convert  --out_orientation RAS --reslice_like %s -rt %s %s %s', alignTo, resample_type, segmentFile, outfile);
eval(str)
    
%%%% read in the nifti and only save all ROIs together
% ni = niftiRead(outfile);
% 
% % map the replacement values
% invals  = [roisToSegment];
% 
% ni.data(~ismember(ni.data,invals))=0;
% ni.data(ismember(ni.data,invals))=1;
% 
% % write out the nifti
% writeFileNifti(ni)

%% make a nifti file from each ROI in aseg.auto.mgz - ROIs to segment

for k=1:length(roisToSegment)
    resample_type= 'weighted';
    
    alignTo = fullfile(dDir,s_info.subj,s_info.anat,[s_info.anatName '.nii']);
    segmentFile = fullfile(dDir,s_info.subj,'freesurfer','mri','aseg.auto.mgz');
    outfile = fullfile(dDir,s_info.subj,'freesurfer','nii',[roisToSegmentNames{k} '.nii.gz']);

    str = sprintf('!mri_convert  --out_orientation RAS --reslice_like %s -rt %s %s %s', alignTo, resample_type, segmentFile, outfile);
    eval(str)

    % read in the nifti
    ni = niftiRead(outfile);

    % map the replacement values
    invals  = [roisToSegment(k)];

    ni.data(~ismember(ni.data,invals))=0;
    ni.data(ismember(ni.data,invals))=1;

    % write out the nifti
    writeFileNifti(ni)
end

%% MAKE ALL FREESURFER LABELS INTO ROIs

%% get ventricle ROI
clear all
close all

dDir = '/biac4/wandell/data/BrainBeat/data';

s = 3;
s_info = bb_subs(s);

% addapted from Franco's s_fe_make_all_freesurfer_labels_into_rois(hemisphere)
% This script creates all the freesurfer lables from parcellation files and 
% creates corresponding nifti and mat ROI compatible with mrDiffusion.
% Copyright by Franco Pestilli Stanford University 2014

% Get the base directory for the data
subjects = {...
%     'yeom_tumor1_5446' ...
    'rosemary2015' ...
    };

fsSubjectsDir = getenv('SUBJECTS_DIR');
annotationFileName = {'aparc','aparc.a2009s'};

clobber = false;

for iSbj = 1:length(subjects)
    fs_subject = subjects{iSbj};
    fsSubjectDir   = fullfile(fsSubjectsDir,fs_subject);
    
    % Create all the necessary label files
    for ia = 1:length(annotationFileName)
        fs_annotationToLabelFiles(fs_subject,annotationFileName{ia},[],fsSubjectsDir);
    end
    
    % File all the label ROIs for this subject
    labelFileNames   = dir(fullfile(fsSubjectDir,'label','*.label'));
    labelRoiName     = cell(length(labelFileNames),1);
    niftiRoiFullPath = cell(length(labelFileNames),1);
    matRoiFullPath  = cell(length(labelFileNames),1);
    for il = 1:length(labelFileNames)
        labelRoiName{il}  = labelFileNames(il).name;
        niftiRoiName      = labelRoiName{il};
        niftiRoiName(niftiRoiName=='.') = '_';
        niftiRoiFullPath{il}  = fullfile(fsSubjectDir,'label',  niftiRoiName);
        matRoiFullPath{il}   = [fullfile(fsSubjectDir,'label',  niftiRoiName),'_smooth3mm_ROI.mat'];
    end
    
    for il = 1:length(labelFileNames)
        if ~(exist([niftiRoiFullPath{il},'_smooth3mm.nii.gz'],'file')==2) || clobber
            fs_labelFileToNiftiRoi(fs_subject,labelRoiName{il},niftiRoiFullPath{il},labelFileNames(il).name(1:2),[],[]);
        else
            fprintf('[%s] Found ROI, skipping: \n%s\n',mfilename,niftiRoiFullPath{il})
        end
        if ~(exist([matRoiFullPath{il}],'file')==2) || clobber
            dtiImportRoiFromNifti([niftiRoiFullPath{il},'_smooth3mm.nii.gz'], matRoiFullPath{il});
        else
            fprintf('[%s] Found ROI, skipping: \n%s\n',mfilename,niftiRoiFullPath{il})
        end
    end
end
