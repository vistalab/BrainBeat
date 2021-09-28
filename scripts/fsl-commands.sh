#!/bin/bash

# Input: $1 - t1w, $2 - fmri 4D data, $3 - Freesurfer T1w, $4 - Freesurfer segmentation nifti, $5 - fmap magnitude, $6 - fmap fieldmap in Hz, $7 - echo spacing in second
# Input nifti should have standard orientation. Use dcm2niix or fslreorient2std to create these niftis

t1w=$1
rest=$2
fs_t1w=$3
fs_aparc=$4
fmap_mag=$5
fmap_frq=$6
esp=$7

# echo $fs_aparc

t1w_outbase="t1w"
rest_outbase="rest"
fs_t1w_outbase="fs_t1w"
fs_aparc_outbase="fs_aparc"
fmap_mag_outbase="fmap_mag"
fmap_frq_outbase="fmap_frq"
xfm_fmap="fmap2epi.txt"
xfm_t1w="t1w2epi.txt"
xfm_fst1w="fst1w2epi.txt"

# unwarp EPI
# convert fieldmap to rad/s and align to the EPI image
$FSLDIR/bin/bet $fmap_mag ${fmap_mag_outbase}_brain -f 0.2
$FSLDIR/bin/bet $rest ${rest_outbase}_brain -f 0.2
$FSLDIR/bin/bet $t1w ${t1w_outbase}_brain -m -f 0.2
$FSLDIR/bin/fslmaths $fmap_frq -mul 6.283 ${fmap_frq_outbase}_rad
$FSLDIR/bin/flirt -in ${fmap_mag_outbase}_brain -ref ${rest_outbase}_brain -dof 6 -cost normmi -omat $xfm_fmap
$FSLDIR/bin/flirt -in ${fmap_frq_outbase}_rad -ref ${rest_outbase}_brain -applyxfm -init $xfm_fmap -out ${fmap_frq_outbase}_rad_resampled
## create brain mask from the anatomical and align to the EPI image
$FSLDIR/bin/flirt -in ${t1w_outbase}_brain -ref ${rest_outbase}_brain -dof 6 -cost normmi -omat $xfm_t1w -out ${t1w_outbase}_brain_resampled
$FSLDIR/bin/flirt -in ${t1w_outbase}_brain_mask -ref ${rest_outbase}_brain -applyxfm -init $xfm_t1w -out ${t1w_outbase}_brain_mask_resampled
## fieldmap correction
$FSLDIR/bin/fugue -i $rest --dwell=$esp --loadfmap=${fmap_frq_outbase}_rad_resampled --mask=${t1w_outbase}_brain_mask_resampled -u ${rest_outbase}_unwarp
#
## register FreeSurfer ROIs. Need to redo alignment because FreeSurfer segmentations are in 1mm iso space instead of the original anatomical space
$FSLDIR/bin/bet ${rest_outbase}_unwarp ${rest_outbase}_unwarp_brain -f 0.2
$FSLDIR/bin/bet $fs_t1w ${fs_t1w_outbase}_brain -f 0.2
$FSLDIR/bin/flirt -in ${fs_t1w_outbase}_brain -ref ${rest_outbase}_unwarp_brain -dof 6 -cost normmi -omat $xfm_fst1w -out ${fs_t1w_outbase}_brain_resampled
$FSLDIR/bin/flirt -in $fs_aparc -ref ${rest_outbase}_unwarp_brain -applyxfm -init $xfm_fst1w -interp nearestneighbour -out ${fs_aparc_outbase}_resampled
