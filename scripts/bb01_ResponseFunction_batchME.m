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
    
    if s==4
        scan_nrs = {[4 5],[6 7]};
    end
    
    for scn = 1%:length(scan_nrs)
        disp(['Calculating response2PPG sub-' int2str(s) ' scan-' int2str(scan_nrs{scn})])

        % first echo
        scan1 = s_info.scan{scan_nrs{scn}(1)};
        scanName1 = s_info.scanName{scan_nrs{scn}(1)};
        fmri1 = fullfile(dDir,subj,scan1,[scanName1 '.nii.gz']);
        niTE1 = niftiRead(fmri1);

        % second echo
        scan2 = s_info.scan{scan_nrs{scn}(2)};
        scanName2 = s_info.scanName{scan_nrs{scn}(2)};
        fmri2 = fullfile(dDir,subj,scan2,[scanName2 '.nii.gz']);
        niTE2 = niftiRead(fmri2);

        % compute the PPG triggered response matrix for the whole brain and save as a nifti:
        [response_matrixS0,t,response_matrix_oddS0,response_matrix_evenS0,response_matrix_stdS0,...
            response_matrixT2s,response_matrix_oddT2s,response_matrix_evenT2s,response_matrix_stdT2s] ...
            = bbResponse2physioME(niTE1,niTE2);

        % safe time T:
        save(fullfile(dDir,subj,scan,[scanName1 '_PPGtrigResponseT']),'t')

        % save average of all heartbeats:
        ni1 = niTE1; ni1.data = response_matrixS0;
        ni1.fname = [scanName1 '_PPGtrigResponse_S0.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
        ni2 = niTE2; ni2.data = response_matrixT2s;
        ni2.fname = [scanName1 '_PPGtrigResponse_T2s.nii.gz'];
        niftiWrite(ni2,fullfile(dDir,subj,scan,ni2.fname))
        clear ni2

        % save standard error of all heartbeats:
        ni1 = niTE1; ni1.data = response_matrix_stdS0;
        ni1.fname = [scanName1 '_PPGtrigResponse_std_S0.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
        ni2 = niTE2; ni2.data = response_matrix_stdT2s;
        ni2.fname = [scanName1 '_PPGtrigResponse_std_T2s.nii.gz'];
        niftiWrite(ni2,fullfile(dDir,subj,scan,ni2.fname))
        clear ni2
        
        % save average of all odd heartbeats:
        ni1 = niTE1; ni1.data = response_matrix_oddS0;
        ni1.fname = [scanName1 '_PPGtrigResponse_odd_S0.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
        ni2 = niTE2; ni2.data = response_matrix_oddT2s;
        ni2.fname = [scanName1 '_PPGtrigResponse_odd_T2s.nii.gz'];
        niftiWrite(ni2,fullfile(dDir,subj,scan,ni2.fname))
        clear ni2
        
        % save average of all even heartbeats:
        ni1 = niTE1; ni1.data = response_matrix_evenS0;
        ni1.fname = [scanName1 '_PPGtrigResponse_even_S0.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1
        ni2 = niTE2; ni2.data = response_matrix_evenT2s;
        ni2.fname = [scanName1 '_PPGtrigResponse_even_T2s.nii.gz'];
        niftiWrite(ni2,fullfile(dDir,subj,scan,ni2.fname))
        clear ni2
    end
    
    for scn = 1%:length(scan_nrs)
        disp(['Calculating PPG locked correlation sub-' int2str(s) ' scan-' int2str(scan_nrs{scn})])
        
        % get scan name
        scan1 = s_info.scan{scan_nrs{scn}(1)};
        scanName1 = s_info.scanName{scan_nrs{scn}(1)};        
        
        % get cod for S0
        ni_odd = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponse_odd_S0.nii.gz']));
        ni_even = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponse_even_S0.nii.gz']));
        
        % compute the cod between the even and odd repetitions of
        % the PPG triggered BOLD signal for the whole brain:
        [out_R_map] = bbCod2physio(ni_odd,ni_even);
        % save as nifti
        ni1 = ni_odd;
        ni1.data = out_R_map;
        ni1.fname = [scanName1 '_S0_codPPG.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1 out_R_map ni_odd ni_even

        % get cod for T2s
        ni_odd = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponse_odd_T2s.nii.gz']));
        ni_even = niftiRead(fullfile(dDir,subj,scan1,[scanName1 '_PPGtrigResponse_even_T2s.nii.gz']));
        
        % compute the cod between the even and odd repetitions of
        % the PPG triggered BOLD signal for the whole brain:
        [out_R_map] = bbCod2physio(ni_odd,ni_even);
        % save as nifti
        ni1 = ni_odd;
        ni1.data = out_R_map;
        ni1.fname = [scanName1 '_T2s_codPPG.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1 out_R_map ni_odd ni_even
        
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

%% Whole brain measurement of reliability
clear all
% close all

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% Select a subject and scan nummer
s_nr = 2;

for scan_nr = 3%[1:3]

    roiNames = {'GM','WM','Ventricles','CSF','Veno'};

    subs = bb_subs(s_nr);
    subj = subs.subj;
    scan = subs.scan{scan_nr};
    scanName = subs.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
%     if ~exist(fmri,'file')
%         clear ni
%         error('filename %s does not exist',fmri)
%     end
%     ni = niftiRead(fmri);
% 
%     % load coregistration matrix for the functionals:
%     load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
%     acpcXform = acpcXform_new; clear acpcXform_new
        
%     % Anatomical
%     anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii.gz']);
%     niAnatomy = niftiRead(anat);

    % COD even/odd
    ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
    ppgR = niftiRead(ppgRname); % correlation with PPG

    % Segmentation file:
    segName = fullfile(dDir,subj,scan,[scanName '_combineSegm.nii.gz']);
    niSeg = niftiRead(segName);

    clear out

    totalNrVox = 0;
    for rr = 1:length(roiNames)
        voxels_inroi = ismember(niSeg.data(:),rr);
        out(rr).roi_ppgR = ppgR.data(voxels_inroi);
        totalNrVox = totalNrVox+length(find(voxels_inroi>0));
    end

    figure('Position',[0 0 700 500])
    subplot(3,2,1),hold on
    r_th = 0.5;
    for rr = 1:length(out)
        bar(rr,100*length(find(out(rr).roi_ppgR>r_th))./length(out(rr).roi_ppgR),'FaceColor',[.5 .5 .5]);
    end
    ylabel(['Percent voxels with R>' num2str(r_th,3)])
    set(gca,'XTick',[1:length(out)],'XTickLabel',roiNames)
    title(['flip angle = ' int2str(subs.scanFA{scan_nr}) ])
    xlim([0 length(out)+1]),ylim([0 100])
    grid on

    % Look at all voxels with R>threshold, and show a pie-chart for how many
    % are in each tissue type.
    subplot(3,2,2),hold on
    nrGray = length(find(ppgR.data(:)>r_th & niSeg.data(:)==1));
    nrWhite = length(find(ppgR.data(:)>r_th & niSeg.data(:)==2));
    nrVentr = length(find(ppgR.data(:)>r_th & niSeg.data(:)==3));
    nrCSF = length(find(ppgR.data(:)>r_th & niSeg.data(:)==4));
    nrVeno = length(find(ppgR.data(:)>r_th & niSeg.data(:)==5));
    
    pie([nrGray,nrWhite,nrVentr,nrCSF,nrVeno],roiNames) 
    axis square
    title(['Location of ' int2str(length(find(ppgR.data(:)>r_th & niSeg.data(:)>0))) ' voxels with R>' num2str(r_th,3)])
    axis off
   
    for kk = 1:5
        subplot(3,5,5+kk)
        hist(ppgR.data(niSeg.data(:)==kk),[0:.1:1])
        h = findobj(gca,'Type','patch');
        h.FaceColor = [.5 .5 .5];
        xlim([0 1])
        title(roiNames{kk})
        xlabel('R'),ylabel('# voxels')
    end
    
    
    % load even and odd responses:
    ni_odd = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_odd.nii.gz']));
    ni_even = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_even.nii.gz']));
    load(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponseT']),'t')

    for kk = 1:5
        subplot(6,5,20+kk),hold on
        
        oddSignals = reshape(ni_odd.data,prod(ni_odd.dim(1:3)),prod(ni_odd.dim(4)));        
        oddThisTissue = oddSignals(ppgR.data(:)>r_th & niSeg.data(:)==kk,:);
        
        plot(t,oddThisTissue)
        xlim([min(t) max(t)])
        ylim([-0.5 .5])
    end
    
    for kk = 1:5
        subplot(6,5,25+kk),hold on
        % rescale:
        oddSignals = ni_odd.data ./ repmat(max(abs(ni_odd.data),[],4),[1,1,1,size(ni_odd.data,4)]);
        % Multiply by correlation size (absolute)
        oddSignals = oddSignals .* abs(repmat(ppgR.data,[1,1,1,size(oddSignals,4)]));
        
        oddSignals = reshape(oddSignals,prod(ni_odd.dim(1:3)),prod(ni_odd.dim(4)));        
        oddThisTissue = oddSignals(ppgR.data(:)>r_th & niSeg.data(:)==kk,:);
        
        plot(t,oddThisTissue,'Color',[.7 .7 .7])
        plot(t,median(oddThisTissue),'k','LineWidth',2)
        
        xlim([min(t) max(t)])
%         ylim([-0.5 .5])
    end
%     set(gcf,'PaperPositionMode','auto')
%     print('-painters','-r300','-dpng',[dDir './figures/reliable/subj' int2str(s_nr) '_scan' int2str(scan_nr)])
%     print('-painters','-r300','-depsc',[dDir './figures/reliable/subj' int2str(s_nr) '_scan' int2str(scan_nr)])
end

