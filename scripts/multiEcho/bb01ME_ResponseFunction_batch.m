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

sub_labels = {'4','5','1'}; 
ses_labels = {'1','1','2'}; 
acq_labels = {'ME','ME','ME'};
run_nrs = {[1],[2],[2]}; % there are a total of 2, 2, 3 runs for these subjects/sessions

for ss = 1%:length(sub_labels) % subjects/ses/acq
    
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    
    for rr = 1:length(run_nrs{ss}) % runs        
        run_nr = run_nrs{ss}(rr);
        
        disp(['Calculating response2PPG sub-' sub_labels{ss} ' run-' int2str(run_nrs{ss}(rr))])
        
        save_dir = fullfile(dDir,'derivatives','brainbeat',['sub-' sub_label],['ses-' ses_label]);
        save_name_base = ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr)];
        
        % first echo
        fmri1_BIDSname = fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'func',...
            ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_echo-1_bold.nii.gz']);
        niTE1 = niftiRead(fmri1_BIDSname);

        % second echo
        fmri2_BIDSname = fullfile(dDir,['sub-' sub_label],['ses-' ses_label],'func',...
            ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_echo-2_bold.nii.gz']);
        niTE2 = niftiRead(fmri2_BIDSname);

        % compute the PPG triggered response matrix for the whole brain and save as a nifti:
        [response_matrixS0,t,response_matrix_oddS0,response_matrix_evenS0,response_matrix_stdS0,...
            response_matrixT2s,response_matrix_oddT2s,response_matrix_evenT2s,response_matrix_stdT2s] ...
            = bbResponse2physioME(niTE1,niTE2);

        % safe time T:
        save(fullfile(save_dir,[save_name_base '_PPGtrigResponseT']),'t')

        % save average of all heartbeats:
        ni1 = niTE1; ni1.data = response_matrixS0;
        ni1.fname = [save_name_base '_PPGtrigResponse_lnS0.nii.gz'];
        niftiWrite(ni1,fullfile(save_dir,ni1.fname))
        clear ni1
        ni2 = niTE2; ni2.data = response_matrixT2s;
        ni2.fname = [save_name_base '_PPGtrigResponse_T2s.nii.gz'];
        niftiWrite(ni2,fullfile(save_dir,ni2.fname))
        clear ni2

        % save standard error of all heartbeats:
        ni1 = niTE1; ni1.data = response_matrix_stdS0;
        ni1.fname = [save_name_base '_PPGtrigResponse_std_lnS0.nii.gz'];
        niftiWrite(ni1,fullfile(save_dir,ni1.fname))
        clear ni1
        ni2 = niTE2; ni2.data = response_matrix_stdT2s;
        ni2.fname = [save_name_base '_PPGtrigResponse_std_T2s.nii.gz'];
        niftiWrite(ni2,fullfile(save_dir,ni2.fname))
        clear ni2
        
        % save average of all odd heartbeats:
        ni1 = niTE1; ni1.data = response_matrix_oddS0;
        ni1.fname = [save_name_base '_PPGtrigResponse_odd_lnS0.nii.gz'];
        niftiWrite(ni1,fullfile(save_dir,ni1.fname))
        clear ni1
        ni2 = niTE2; ni2.data = response_matrix_oddT2s;
        ni2.fname = [save_name_base '_PPGtrigResponse_odd_T2s.nii.gz'];
        niftiWrite(ni2,fullfile(save_dir,ni2.fname))
        clear ni2
        
        % save average of all even heartbeats:
        ni1 = niTE1; ni1.data = response_matrix_evenS0;
        ni1.fname = [save_name_base '_PPGtrigResponse_even_lnS0.nii.gz'];
        niftiWrite(ni1,fullfile(save_dir,ni1.fname))
        clear ni1
        ni2 = niTE2; ni2.data = response_matrix_evenT2s;
        ni2.fname = [save_name_base '_PPGtrigResponse_even_T2s.nii.gz'];
        niftiWrite(ni2,fullfile(save_dir,ni2.fname))
        clear ni2
    end
    
%     for scn = 1:length(scan_nrs)
%         disp(['Calculating PPG locked correlation sub-' int2str(s) ' scan-' int2str(scan_nrs{scn})])
%         
%         % get scan name
%         scan1 = s_info.scan{scan_nrs{scn}(1)};
%         scanName1 = s_info.scanName{scan_nrs{scn}(1)};        
%         
%         % get cod for S0
%         ni_odd = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponse_odd_S0.nii.gz']));
%         ni_even = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponse_even_S0.nii.gz']));
%         
%         % compute the cod between the even and odd repetitions of
%         % the PPG triggered BOLD signal for the whole brain:
%         [out_R_map] = bbCod2physio(ni_odd,ni_even);
%         % save as nifti
%         ni1 = ni_odd;
%         ni1.data = out_R_map;
%         ni1.fname = [scanName1 '_S0_codPPG.nii.gz'];
%         niftiWrite(ni1,fullfile(dDir,subj,scan1,ni1.fname))
%         clear ni1 out_R_map ni_odd ni_even
% 
%         % get cod for T2s
%         ni_odd = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponse_odd_T2s.nii.gz']));
%         ni_even = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponse_even_T2s.nii.gz']));
%         
%         % compute the cod between the even and odd repetitions of
%         % the PPG triggered BOLD signal for the whole brain:
%         [out_R_map] = bbCod2physio(ni_odd,ni_even);
%         % save as nifti
%         ni1 = ni_odd;
%         ni1.data = out_R_map;
%         ni1.fname = [scanName1 '_T2s_codPPG.nii.gz'];
%         niftiWrite(ni1,fullfile(dDir,subj,scan1,ni1.fname))
%         clear ni1 out_R_map ni_odd ni_even
%         
%     end
end

%%
figure('Position',[0 0 300 400])
subplot(3,1,1),hold on
plot(t,squeeze(response_matrixS0(27,50,20,:)))
subplot(3,1,2),hold on
plot(t,-1./squeeze(response_matrixT2s(27,50,20,:)))
subplot(3,1,3),hold on
plot(t,squeeze(response_matrixT2s(27,50,20,:)))
figName = fullfile(bbPath,'local','ME_ExampleAvTrace');

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',figName)
print('-painters','-r300','-depsc',figName)
