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
    
    for scan_nr = [1:9]%1:length(s_info.scan)
        scan = s_info.scan{scan_nr};
        scanName = s_info.scanName{scan_nr};

        ni_odd = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_odd.nii.gz']));
        ni_even = niftiRead(fullfile(dDir,subj,scan,[scanName '_PPGtrigResponse_even.nii.gz']));
        
        % compute the correlation between the even and odd repetitions of
        % the PPG triggered BOLD signal for the whole brain:
        [out_r_map] = bbCorrelate2physio(ni_odd,ni_even);
        % save as nifti
        ni1 = ni_odd;
        ni1.data = out_r_map;
        ni1.fname = [scanName '_corrPPG.nii.gz'];
        niftiWrite(ni1,fullfile(dDir,subj,scan,ni1.fname))
        clear ni1 out_r_map

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

%% Whole brain measurement of reliability

% Select a subject and scan nummer
s_nr = 4;

for scan_nr = [1:9]

    roisToSegment = {[...
        4 % left lateral ventricle
        5 % left inferior lateral ventricle
        14 % 3rd ventricle
        15 % 4th ventricle
        24 % CSF
        31 % left choroid plexus
        43 % right lateral ventricle
        44 % right inferior lateral ventricle
        63 % right choroid plexus
        72],... % 5th ventricle
        [2 % left white matter
        41],... % right white matter
        [3 % left gray matter
        42]}; % right gray matter
    roiNames = {'CSF','WM','GM','Veno'};

    subs = bb_subs(s_nr);
    subj = subs.subj;
    scan = subs.scan{scan_nr};
    scanName = subs.scanName{scan_nr};
    fmri = fullfile(dDir,subj,scan,[scanName '.nii.gz']);
    if ~exist(fmri,'file')
        clear ni
        error('filename %s does not exist',fmri)
    end
    ni = niftiRead(fmri);

    % load coregistration matrix for the functionals:
    load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
    acpcXform = acpcXform_new; clear acpcXform_new

    % Get the MRVenogram:
    niVeno = niftiRead(fullfile(dDir,subj,subs.veno,[subs.venoName '.nii.gz']));
    % load coregistration matrix (for the venogram):
    xf_veno=load(fullfile(dDir,subj,subs.veno,[subs.venoName 'AcpcXform.mat']));

    % Anatomical
    anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii.gz']);
    niAnatomy = niftiRead(anat);

    % COD even/odd
    ppgRname = fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz']);
    ppgR = niftiRead(ppgRname); % correlation with PPG

    % Segmentation file:
    segName = fullfile(dDir,subj,scan,[scanName '_r_aseg_auto.nii.gz']);
    niSeg = niftiRead(segName);

    % Veno file:
    venoName = fullfile(dDir,subj,scan,[scanName '_r_veno.nii.gz']);
    niVeno = niftiRead(venoName);

    % % Get voxels with mean signal above brain_th
    % brain_th = 50;
    % brain_vect = mean(ni.data(:,:,:,4:end),4);
    % brain_vect = brain_vect>brain_th;

    clear out

    totalNrVox = 0;
    for rr = 1:length(roisToSegment)
        voxels_inroi = ismember(niSeg.data(:),roisToSegment{rr});
        out(rr).roi_ppgR = ppgR.data(voxels_inroi);
        totalNrVox = totalNrVox+length(find(voxels_inroi>0));
    end

    voxels_inroi = niVeno.data(:)>1000;
    out(rr+1).roi_ppgR = ppgR.data(voxels_inroi);

    r_th = .3;

%     figure('Position',[0 0 180 200])
%     hold on
%     for rr = 1:length(out)
%         bar(rr,100*length(find(out(rr).roi_ppgR>r_th))./length(out(rr).roi_ppgR));
%     end
%     ylabel(['Percent voxels with R>' num2str(r_th,3)])
%     set(gca,'XTick',[1:length(out)],'XTickLabel',roiNames)
%     title(['flip angle = ' int2str(subs.scanFA{scan_nr}) ])
%     xlim([0 length(out)+1]),ylim([0 100])
%     grid on
%     
%     set(gcf,'PaperPositionMode','auto')
%     print('-painters','-r300','-dpng',[dDir './figures/reliable/subj' int2str(s_nr) '_scan' int2str(scan_nr)])
%     print('-painters','-r300','-depsc',[dDir './figures/reliable/subj' int2str(s_nr) '_scan' int2str(scan_nr)])

    % Look at all voxels with R>threshold, and show a pie-chart for how many
    % are in each tissue type.
    r_th = .3;
    veno_th = 1000;
    func_mean = mean(ni.data(:,:,:,5:end),4);
    brain_th = mean(func_mean(:));
    
    allR = ppgR.data(:)>r_th & func_mean(:)>brain_th;
    allVeno = allR(:) & niVeno.data(:)>veno_th & func_mean(:)>brain_th;
    nrCSF = length(find(allR & ismember(niSeg.data(:),roisToSegment{1}) & ~allVeno));
    nrWhite = length(find(allR & ismember(niSeg.data(:),roisToSegment{2}) & ~allVeno));
    nrGray = length(find(allR & ismember(niSeg.data(:),roisToSegment{3}) & ~allVeno));
    nrVeno = length(find(allVeno>0));
    nrOther = length(find(allR>0)) - sum([nrCSF,nrWhite,nrGray,nrVeno]);

    figure('Position',[0 0 200 200])
    labels = {'Other','CSF','White','Gray','Veno'};
    pie([nrOther,nrCSF,nrWhite,nrGray,nrVeno],labels) 
    title(['flip angle = ' int2str(subs.scanFA{scan_nr})])

    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',[dDir './figures/reliable/pie_subj' int2str(s_nr) '_scan' int2str(scan_nr)])
    print('-painters','-r300','-depsc',[dDir './figures/reliable/pie_subj' int2str(s_nr) '_scan' int2str(scan_nr)])

end

