BrainBeat
=========



## Notes on BIDS formatting

### Slice timing
The _bold.json file contains a field with the SliceTiming. This is essential to get the timing accuracy for every slice.

### Physio data
Physio data include PPG, RESP in all subjects and cardiac ECG measurements in a few subjects. Sampling frequencies differ across subjects and these data have to be stored in different .tsv.gz files in the BIDS format.
PPG physio files have to have the same name as the recording plus a _recording-PPG label.
RESP physio files have to have the same name as the recording plus a _recording-RESP label.
ECG physio files have to have the same name as the recording plus a _recording-ECG label.

## Dependencies

- vistasoft: https://github.com/vistalab/vistasoft
- bids-matlab is located in utilities/external for reading tsv and json files (source: https://github.com/bids-standard/bids-matlab)
- spm12: https://github.com/spm/spm12 (addpath to this and run spm('Defaults','fmri'))

Matlab toolboxes:
- signal processing toolbox
- statistics toolbox
- image processing toolbox

## Contact
Dora Hermes: 
Brian Wandell:

