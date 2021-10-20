clear all
close all

%% s_bbIntro
%
% Script to explain how to open and do preliminary brain beat analyses


%% Base data directory 

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';


%% Preprocess the data: PPG time series and correlation

sub_labels = {'1'}; 
ses_labels = {'2'}; 
acq_labels = {'4mmFA48'};
run_nrs = {[2]};

for ss = 1:length(sub_labels) % subjects/ses/acq
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    
    for rr = 1:length(run_nrs{ss}) % runs
        
        run_nr = run_nrs{ss}(rr);
        
        fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
            ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
        fmri_name = fullfile(dDir,fmri_BIDSname);
        save_dir = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label]);
        if ~exist(save_dir,'dir')
            disp(['creating output directory ' save_dir])
            mkdir(save_dir)
        end
        save_name_base = ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)];

        if ~exist(fmri_name,'file')
            clear ni
            error('filename %s does not exist',fmri_name)
        end
        ni = niftiRead(fmri_name);

        disp(['Calculating PPG triggered responses for ' fmri_BIDSname])

        % compute the PPG triggered response matrix for the whole brain and save as a nifti:
        [response_matrix,t,response_matrix_odd,response_matrix_even,response_matrix_std] ...
            = bbResponse2physio(ni);

        % safe time T:
        save(fullfile(save_dir,[save_name_base '_PPGtrigResponseT']),'t')

        % save average of all heartbeats:
        ni1 = ni;
        ni1.data = response_matrix;
        ni1.fname = [save_name_base '_PPGtrigResponse.nii.gz'];
        niftiWrite(ni1,fullfile(save_dir,ni1.fname))
        clear ni1

        % save standard deviation of all heartbeats:
        ni1 = ni;
        ni1.data = response_matrix_std;
        ni1.fname = [save_name_base '_PPGtrigResponse_std.nii.gz'];
        niftiWrite(ni1,fullfile(save_dir,ni1.fname))
        clear ni1
        
        % save average of all odd heartbeats:
        ni1 = ni;
        ni1.data = response_matrix_odd;
        ni1.fname = [save_name_base '_PPGtrigResponse_odd.nii.gz'];
        niftiWrite(ni1,fullfile(save_dir,ni1.fname))
        clear ni1
        
        % save average of all even heartbeats:
        ni1 = ni;
        ni1.data = response_matrix_even;
        ni1.fname = [save_name_base '_PPGtrigResponse_even.nii.gz'];
        niftiWrite(ni1,fullfile(save_dir,ni1.fname))
        clear ni1
    end
    

    for rr = 1:length(run_nrs{ss}) % runs
        
        run_nr = run_nrs{ss}(rr);
        
        disp(['Calculating PPG locked reliability (COD) for ' fmri_BIDSname])
        
        fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
            ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
        fmri_name = fullfile(dDir,fmri_BIDSname);
        save_dir = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label]);
        save_name_base = ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)];
        
        ni_odd = niftiRead(fullfile(save_dir,[save_name_base '_PPGtrigResponse_odd.nii.gz']));
        ni_even = niftiRead(fullfile(save_dir,[save_name_base '_PPGtrigResponse_even.nii.gz']));
        
        % compute the COD between the even and odd repetitions of
        % the PPG triggered BOLD signal for the whole brain:
        [out_R_map] = bbCod2physio(ni_odd,ni_even);
        
        % save as nifti
        ni2 = ni_odd;
        ni2.data = out_R_map;
        ni2.fname = [save_name_base '_codPPG.nii.gz'];
        niftiWrite(ni2,fullfile(save_dir,ni2.fname))
        clear ni2 out_R_map
       
    end
end



%%
%% This has not been adjusted to BIDS yet...
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

