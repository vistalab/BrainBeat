function data_colors_rgb = bbData2Colors(data_in)
% 
% data_in should scale [-1:1]
%
% % example data and figure
% data_in = [0 0;0 1;1 0;1 1;0 -1;-1 0;-1 -1;-1 1;1 -1];
% data_in = [data_in;.2 * data_in; .4 * data_in; .6 * data_in; .8 * data_in];
% data_colors_rgb = bbData2Colors(data_in);
% figure,hold on
% for k=1:length(data_in)
%     plot(data_in(k,1),data_in(k,2),'k.','MarkerSize',40)
%     plot(data_in(k,1),data_in(k,2),'.','Color',data_colors_rgb(k,:),'MarkerSize',35)
% end
% print('-painters','-r300','-dpng',['./figures/test_pca/colorspace01'])

% make complex:
data = complex(data_in(:,1),data_in(:,2));

% convert complex data to hsv colors, can be done in different ways:
% rotate the hue
rot_f =  0;%pi/4; % zero - no rotation

% data_colors_hsv = [(pi+angle(data)+rot_f)/(2*pi) abs(data)/sqrt(2) abs(data)/sqrt(2)];
data_colors_hsv = [mod(angle(data)+pi+rot_f,2*pi)/(2*pi) max(abs(data_in),[],2) max(abs(data_in),[],2)];
% data_colors_hsv = [mod(angle(data)+pi+rot_f,2*pi)/(2*pi) abs(data)/sqrt(2) abs(data)/sqrt(2)];

% data_colors_hsv = [(pi+angle(data))/(2*pi) abs(data)/sqrt(2) .5*ones(size(data))];
% data_colors_hsv = [(pi+angle(data))/(2*pi) abs(data)/sqrt(2) .9*ones(size(data))];

% get rgb values instead of hsv
data_colors_rgb = hsv2rgb(data_colors_hsv);




