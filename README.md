BrainBeat
=========

This repository contains functions to do the analyses and generate figures from the manuscript titled:

**Measuring brain beats: cardiac-aligned fast fMRI signals.** <br>
Dora Hermes, Hua Wu, Adam Kerr, Brian Wandell, <br> 
*Human Brain Mapping 2022* https://doi.org/10.1002/hbm.26128

Please cite this work when using the code.

**Abstract**  Blood and cerebrospinal fluid (CSF) pulse and flow throughout the brain, driven by the cardiac cycle. These fluid dynamics, which are essential to healthy brain function, are characterized by several noninvasive magnetic resonance imaging (MRI) methods. Recent developments in fast MRI, specifically simultaneous multislice acquisition methods, provide a new opportunity to rapidly and broadly assess cardiac-driven flow, including CSF spaces, surface vessels and parenchymal vessels. We use these techniques to assess blood and CSF flow dynamics in brief (3.5 min) scans on a conventional 3 T MRI scanner in five subjects. Cardiac pulses are measured with a photoplethysmography (PPG) on the index finger, along with functional MRI (fMRI) signals in the brain. We, retrospectively, align the fMRI signals to the heartbeat. Highly reliable cardiac-gated fMRI temporal signals are observed in CSF and blood on the timescale of one heartbeat (test–retest reliability within subjects R2 > 50%). In blood vessels, a local minimum is observed following systole. In CSF spaces, the ventricles and subarachnoid spaces have a local maximum following systole instead. Slower resting-state scans with slice timing, retrospectively, aligned to the cardiac pulse, reveal similar cardiac-gated responses. The cardiac-gated measurements estimate the amplitude and phase of fMRI pulsations in the CSF relative to those in the arteries, an estimate of the local intracranial impedance. Cardiac aligned fMRI signals can provide new insights about fluid dynamics or diagnostics for diseases where these dynamics are important.

## Download data
Data related to this manuscript are located on [OpenNeuro](https://openneuro.org/datasets/ds004213).

Download data and save them in the `/local/BrainBeat/` folder in this repository.

## Recreate figures from the manuscript
Add appropriate paths (see notes on dependencies below)
Run: `main_createFigures.m`

Figures will be saved in: `/BrainBeat/local/BrainBeat/derivatives/figures/`

Here is an example:  ![image](https://github.com/user-attachments/assets/bb325bc4-a4f7-4994-9846-8827deb6b219)

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

