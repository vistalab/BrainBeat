BrainBeat
=========

This repository contains functions to do the analyses and generate figures from the manuscript titled:

**Measuring brain beats: cardiac-aligned fast fMRI signals.**
Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell, *Human Brain Mapping 2022*

Please cite this work when using the code.

## Download data
Data related to this manuscript are located on [OpenNeuro](https://openneuro.org/datasets/ds004213).

Download data and save them in the `/local/BrainBeat/` folder in this repository.

## Recreate figures from the manuscript
Add appropriate paths (see notes on dependencies below)
Run: `main_createFigures.m`

Figures will be saved in: `/BrainBeat/local/BrainBeat/derivatives/figures/`

## Dependencies

General toolboxes
- [vistasoft](https://github.com/vistalab/vistasoft)
- [bids-matlab](https://github.com/bids-standard/bids-matlab)
- [spm12](https://github.com/spm/spm12), never use genpath when adding SPM in your matlab path, always addpath to the main SMP folder and run spm('Defaults','fmri'))

Matlab toolboxes:
- signal processing toolbox
- statistics toolbox
- image processing toolbox

## Notes on data and BIDS formatting

Anatomical MRI images were de-identified using the [mri_reface](https://www.nitrc.org/projects/mri_reface) function from Christopher Schwarz.

### Slice timing
The _bold.json file contains a field with the SliceTiming. This is essential to get the timing accuracy for every slice.

### Physio data
Physio data include PPG, RESP in all subjects and cardiac ECG measurements in a few subjects. Sampling frequencies differ across subjects and these data have to be stored in different .tsv.gz files in the BIDS format.
- PPG physio files have to have the same name as the recording plus a _recording-PPG label.
- RESP physio files have to have the same name as the recording plus a _recording-RESP label.
- ECG physio files have to have the same name as the recording plus a _recording-ECG label.


## Contact
Dora Hermes: dorahermes@gmail.com

