
% CSF and blood map
cm_csf = customcolormap([0 .45 .55 1], [.7 1 .6;.5 1 .9;.6 .1 .6;1 .9 .8], 200);
cm_blood = customcolormap([0 .45 .55 1], [.4 1 1;.4 .4 1; 1 .4 .4;1 1 .4], 200);
% add gray
gray_vect = .2*ones(200,3);
cm2D_csf = zeros(100,size(cm_csf,1),3);
cm2D_blood = zeros(100,size(cm_blood,1),3);
for kk = 1:100
    cm2D_csf(kk,:,:) = cm_csf*kk/100 + gray_vect*(100-kk)/100;
    cm2D_blood(kk,:,:) = cm_blood*kk/100 + gray_vect*(100-kk)/100;
end

figure,hold on
for kk = 1:size(cm2D_csf,1)
    for mm = 1:size(cm2D_csf,2)
        plot(mm,kk,'.','MarkerSize',25,'Color',cm2D_csf(kk,mm,:))
    end
end
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','colormap_csf'))

figure,hold on
for kk = 1:size(cm2D_blood,1)
    for mm = 1:size(cm2D_blood,2)
        plot(mm,kk,'.','MarkerSize',25,'Color',cm2D_blood(kk,mm,:))
    end
end
set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',fullfile(dDir,'derivatives','brainbeat','group','colormap_blood'))
