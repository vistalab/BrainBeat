
%% plot Resp

all_phys = [];

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

for ss = 1:length(sub_labels)    
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    rr = 1;    
    run_nr = run_nrs{ss}(rr);

    fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
    fmri_name = fullfile(dDir,fmri_BIDSname);
    ni = niftiRead(fmri_name);

    % get physio stuff we need:
    all_phys(ss).physio = physioCreate('nifti',ni);
    all_phys(ss).resp_onsets = physioGet(all_phys(ss).physio,'resp peaks');
end


%%
figure('Position',[0 0 800 800])

for ss = 1:length(sub_labels) 

    subplot(6,2,ss*2-1),hold on

    resp_srate = all_phys(ss).physio.resp.srate;
    resp_plot = zscore(all_phys(ss).physio.resp.data);
    plot([1:length(all_phys(ss).physio.resp.data)]/resp_srate,resp_plot,'Color',[.2 .2 .2])
    plot(all_phys(ss).resp_onsets,resp_plot(round(all_phys(ss).resp_onsets*resp_srate)),'k*')
    xlim([0 all_phys(ss).resp_onsets(end)])
%     xlim([79 92])
    xlabel('time (s)')

    subplot(6,2,ss*2),hold on
    plot(all_phys(ss).resp_onsets(2:end),diff(all_phys(ss).resp_onsets),'Color',[.5 .5 .5])
    xlabel('time (s)')
    xlim([0 all_phys(ss).resp_onsets(end)])

end

subplot(6,2,5),
ylabel('respiratory signal (z-score)')

subplot(6,2,6),
ylabel('breath-breath interval (s)')

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','RESP_allsub'))
print('-painters','-r300','-depsc2',fullfile(dDir,'derivatives','brainbeat','RESP_allsub'))


%% get Heart Rate Variability

all_phys = [];

sub_labels = {'1','2','3','4','5','1'}; 
ses_labels = {'1','1','1','1','1','2'}; 
acq_labels = {'4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48','4mmFA48'};
run_nrs = {[1],[1],[1],[1],[1],[1]};

for ss = 1:length(sub_labels)    
    sub_label = sub_labels{ss};
    ses_label = ses_labels{ss};
    acq_label = acq_labels{ss};
    rr = 1;    
    run_nr = run_nrs{ss}(rr);

    fmri_BIDSname = fullfile(['sub-' sub_label],['ses-' ses_label],'func',...
        ['sub-' sub_label '_ses-' ses_label '_acq-' acq_label '_run-' int2str(run_nr) '_bold.nii.gz']);
    fmri_name = fullfile(dDir,fmri_BIDSname);
    ni = niftiRead(fmri_name);

    % get physio stuff we need:
    all_phys(ss).physio = physioCreate('nifti',ni);
    all_phys(ss).ppg_onsets = physioGet(all_phys(ss).physio,'ppg peaks');
end

%% plot HRV

cl_use = lines(6);

figure

% just subsequent differences (RRs)
for ss = 1:6
    subplot(2,4,1:3),hold on
    plot(all_phys(ss).ppg_onsets(2:end),diff(all_phys(ss).ppg_onsets),'Color',cl_use(ss,:));
    
    subplot(2,4,4),hold on    
    b = bar(ss,mean(diff(all_phys(ss).ppg_onsets)));
    set(b,'FaceColor',cl_use(ss,:))
    eb = errorbar(ss,mean(diff(all_phys(ss).ppg_onsets)),std(diff(all_phys(ss).ppg_onsets)));
    set(eb,'Color','k')   
end

subplot(2,4,1:3),hold on
ylabel('P-P interval (s)')
xlabel('Time (s)')
xlim([0 max(all_phys(1).ppg_onsets)])

subplot(2,4,4),hold on
set(gca,'XTick',[1:6])
xlabel('Subject #'),ylabel('P-P interval (s)')

% RMSSD
int_length = 20; % seconds to calculate RMSSD
all_hrv = zeros(6,1);
for ss = 1:6
    subplot(2,4,5:7),hold on
    these_int = 1:int_length:max(all_phys(ss).ppg_onsets);
    this_hrv = zeros(length(these_int),1);
    for kk = 1:length(these_int)
        % get 10 seconds
        these_beats = all_phys(ss).ppg_onsets(all_phys(ss).ppg_onsets>these_int(kk) & all_phys(ss).ppg_onsets<=these_int(kk)+int_length);
        % get RR (peak to peak difference),
        RR = diff(these_beats);
        % difference in ms between RRs
        these_diff = diff(RR)*1000; 
        this_hrv(kk) = sqrt(mean(these_diff.^2));
    end
    plot(these_int-1+int_length,this_hrv,'Color',cl_use(ss,:))
    all_hrv(ss) = mean(this_hrv);
    
    subplot(2,4,8),hold on
    b = bar(ss,all_hrv(ss));
    set(b,'FaceColor',cl_use(ss,:))   
    eb = errorbar(ss,all_hrv(ss),std(this_hrv));
    set(eb,'Color','k')   
end

subplot(2,4,5:7),hold on
title('Root Mean Square of Successive Differences between heartbeats')
xlabel(['Time (s)']),ylabel('RMSSD (ms)')
xlim([0 max(all_phys(1).ppg_onsets)])
% b.CData = cl_use;
subplot(2,4,8),hold on
set(gca,'XTick',[1:6])
ylabel('RMSSD (ms)')
xlabel('Subject #')

print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','HRVinPPG'))
print('-painters','-r300','-depsc2',fullfile(dDir,'derivatives','brainbeat','group','HRVinPPG'))
