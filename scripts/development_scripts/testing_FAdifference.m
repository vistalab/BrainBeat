ni_25 = niftiRead([dDir '/sourcedata/20141017_1242/6_1_mux8fov4_r1_25s_4mmFA25/8202_6_1.nii.gz']);


ni_34 = niftiRead([dDir '/sourcedata/20141017_1242/5_1_mux8fov4_r1_25s_4mm/8202_5_1.nii.gz']);


ni_48 = niftiRead([dDir '/sourcedata/20141017_1242/7_1_mux8fov4_r1_25s_4mmFA48/8202_7_1.nii.gz']);


%%


mean_25 = mean(ni_25.data(:,:,:,5:end),4);
mean_34 = mean(ni_34.data(:,:,:,5:end),4);
mean_48 = mean(ni_48.data(:,:,:,5:end),4);


%%
aa = mean_34./mean_25;
bb = mean_48./mean_25;
cc = mean_48./mean_34;

figure

subplot(1,3,1),imagesc(aa(:,:,22)',[0 2])
subplot(1,3,2),imagesc(bb(:,:,22)',[0 2])
subplot(1,3,3),imagesc(cc(:,:,22)',[0 2])

%%
figure
subplot(1,3,1),imagesc(squeeze(aa(27,:,:))',[0 2])
axis xy
subplot(1,3,2),imagesc(squeeze(bb(27,:,:))',[0 2])
axis xy
subplot(1,3,3),imagesc(squeeze(cc(27,:,:))',[0 2])
axis xy

%%
acpcXform_new = load([dDir '20141017_1242/5_1_mux8fov4_r1_25s_4mm/8202_5_1AcpcXform_new.mat']);
acpcXform = acpcXform_new.acpcXform_new;

niAnatomy = niftiRead([dDir 'sourcedata/20141017_1242/9_1_T1w_1mm_sag/8202_9_1.nii.gz']);

thisPlot = ni_25;
thisPlot.data = bb-1;

%% Subject 1 Sagittal slices T1 + COD
curPos = [-10,1,-25]; 
sliceThisDim = 3; 
imDims = [-90 -120 -120; 90 130 90];
overlayPlot = thisPlot;
cod_th = 2;

for kk = -5
    curPos(1) = kk;
%     bbOverlayFuncAnat(overlayPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,2.5,1,.6);
    bbOverlayFuncAnat(overlayPlot,niAnatomy,acpcXform,sliceThisDim,imDims,curPos,.5,1,-1);
    title(['slice ' int2str(kk) ', R^2>' num2str(cod_th,3)])
    set(gcf,'PaperPositionMode','auto')    
%     print('-painters','-r300','-dpng',[dDir '/figures/reliable/sub-' int2str(s_nr) '_scan-' int2str(scan_nr) '_orient' int2str(sliceThisDim) '_slice' int2str(kk)])
end
