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

for s = 5
    s_info = bb_subs(s);
    subj = s_info.subj;
    for scan_nr = [2]%:length(s_info.scan)
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
    
    for scan_nr = [2]%1:length(s_info.scan)
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


%%
%% Mean signals in different tissue types
%%

clear all
close all

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% Select a subject and scan nummer
s_nr = 2;

figure('Position',[0 0 1000 400])
scan_nrs = [3];
for ss = 1:length(scan_nrs)
    
    scan_nr = scan_nrs(ss);
    
    roiNames = {'G','W','V','C','Veno'};

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
 
%     % load coregistration matrix for the functionals:
%     load(fullfile(dDir,subj,scan,[scanName 'AcpcXform_new.mat']))
%     acpcXform = acpcXform_new; clear acpcXform_new
        
%     % Anatomical
%     anat      = fullfile(dDir,subj,subs.anat,[subs.anatName '.nii.gz']);
%     niAnatomy = niftiRead(anat);

    % Segmentation file:
    segName = fullfile(dDir,subj,scan,[scanName '_combineSegm.nii.gz']);
    niSeg = niftiRead(segName);

    clear out

    totalNrVox = 0;
    for rr = 1:length(roiNames)
        voxels_inroi = ismember(niSeg.data(:),rr);
        niMean = mean(ni.data(:,:,:,10:end));
        out(rr).mean = niMean(voxels_inroi);
    end

    
    subplot(1,length(scan_nrs),ss),hold on
    r_th = 0.3;
    for rr = 1:length(out)
        bar(rr,mean(out(rr).mean),'FaceColor',[.5 .5 .5]);
        errorbar(rr,mean(out(rr).mean),2*std(out(rr).mean)./sqrt(length(out(rr).mean)),'k');
    end
    
    set(gca,'XTick',[1:length(out)],'XTickLabel',roiNames,'YTick',[1500:500:4000],'YTickLabel',[])
    title(['flip angle = ' int2str(subs.scanFA{scan_nr}) ])
    xlim([0 length(out)+1])%,ylim([1500 4000])
    grid on
    
    
end

subplot(1,length(scan_nrs),1),hold on
ylabel(['Mean signal'])
set(gca,'YTick',[1500:500:4000],'YTickLabel',{'','2000','','3000','','4000'})
set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',[dDir './figures/reliable/MeanSig_subj' int2str(s_nr)])
% print('-painters','-r300','-depsc',[dDir './figures/reliable/MeanSig_subj' int2str(s_nr)])


%%
%% Reliability in SPM segmentation
clear all
% close all

dDir = '/Volumes/DoraBigDrive/data/BrainBeat/data/';

% Select a subject and scan nummer
s_nr = 5;

for scan_nr = 1%[1:3]

    roiNames = {'GM','WM','CSF'};
    roiColormap = ([.5 .5 .5;1 1 1;0 .4 .8]);

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

    % COD even/odd -
    ppgR = niftiRead(fullfile(dDir,subj,scan,[scanName '_codPPG.nii.gz'])); 

    % SPM segmentation
    niSPM = niftiRead(fullfile(dDir,subj,scan,[scanName '_spmSeg.nii.gz']));

    clear out

    totalNrVox = 0;
    for rr = 1:length(roiNames)
        voxels_inroi = ismember(niSPM.data(:),rr);
        out(rr).roi_ppgR = ppgR.data(voxels_inroi);
        totalNrVox = totalNrVox+length(find(voxels_inroi>0));
    end

    figure('Position',[0 0 700 500])
    subplot(3,5,1),hold on
    r_th = 0.5;
    for rr = 1:length(out)
        bar(rr,100*length(find(out(rr).roi_ppgR>r_th))./length(out(rr).roi_ppgR),'FaceColor',roiColormap(rr,:));
    end
    ylabel(['Percent voxels with R>' num2str(r_th,3)])
    set(gca,'XTick',[1:length(out)],'XTickLabel',roiNames)
    title(['flip angle = ' int2str(subs.scanFA{scan_nr}) ])
    xlim([0 length(out)+1]),ylim([0 100])
    grid on

    % Look at all voxels with R>threshold, and show a pie-chart for how many
    % are in each tissue type.
    subplot(3,2,2),hold on
    nrGray = length(find(ppgR.data(:)>r_th & niSPM.data(:)==1));
    nrWhite = length(find(ppgR.data(:)>r_th & niSPM.data(:)==2));
    nrCSF = length(find(ppgR.data(:)>r_th & niSPM.data(:)==3));
    
    pie([nrGray,nrWhite,nrCSF],roiNames) 
    axis square
    colormap(roiColormap)
    title(['Location of ' int2str(length(find(ppgR.data(:)>r_th & niSPM.data(:)>0))) ' voxels with R>' num2str(r_th,3)])
    axis off
   
    for kk = 1:length(out)
        subplot(3,length(out),length(out)+kk)
        hist(ppgR.data(niSPM.data(:)==kk),[0:.1:1])
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

    for kk = 1:length(out)
        subplot(3,length(out),2*length(out)+kk),hold on
        
        oddSignals = reshape(ni_odd.data,prod(ni_odd.dim(1:3)),prod(ni_odd.dim(4)));        
        oddThisTissue = oddSignals(ppgR.data(:)>r_th & niSPM.data(:)==kk,:);
        
        plot(t,100*oddThisTissue)
        xlim([min(t) max(t)])
        ylabel('percent change')
    end
    
    set(gcf,'PaperPositionMode','auto')
    print('-painters','-r300','-dpng',[dDir './figures/reliable/subj' int2str(s_nr) '_scan' int2str(scan_nr)])
    print('-painters','-r300','-depsc',[dDir './figures/reliable/subj' int2str(s_nr) '_scan' int2str(scan_nr)])
end