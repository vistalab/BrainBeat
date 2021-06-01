BrainBeat
=========



### Notes on BIDS formatting

Physio data include PPG, RESP in all subjects and cardiac measurements in a few subjects. Sampling frequencies differ across subjects and these data have to be stored in different .tsv.gz files. 
PPG physio files have to have the same name as the recording plus a _recording-PPG label.
RESP physio files have to have the same name as the recording plus a _recording-RESP label.


### Dependencies

bids-matlab is located in the external folder for reading tsv and json files.

Matlab toolboxes:
- signal processing toolbox


### Contact
