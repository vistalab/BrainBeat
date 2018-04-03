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
    
    for scan_nr = [2 4 5 7 8]%1:length(s_info.scan)
        scan = s_info.scan{scan_nr};
        scanName = s_info.scanName{scan_nr};

        ni_odd = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_odd.nii.gz']));
        ni_even = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_even.nii.gz']));
        
%         % compute the correlation between the even and odd repetitions of
%         % the PPG triggered BOLD signal for the whole brain:
%         [out_r_map] = bbCorrelate2physio(ni_odd,ni_even);
%         % save as nifti
%         ni1 = ni_odd;
%         ni1.data = out_r_map;
%         ni1.fname = [scanName '_corrPPG.nii.gz'];
%         niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
%         clear ni1 out_r_map

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
