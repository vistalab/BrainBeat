data_in = [0 0;0 1;1 0;1 1;0 -1;-1 0;-1 -1;-1 1;1 -1];
data_in = [data_in;.5 * data_in];

% make complex:
data = complex(data_in(:,1),data_in(:,2));

data_colors_hsv = [(pi+angle(data))/(2*pi) abs(data)/sqrt(2) abs(data)/sqrt(2)];
% data_colors_hsv = [(pi+angle(data))/(2*pi) abs(data)/sqrt(2) .5*ones(size(data))];
% data_colors_hsv = [(pi+angle(data))/(2*pi) abs(data)/sqrt(2) .9*ones(size(data))];
data_colors_rgb = hsv2rgb(data_colors_hsv);

figure,hold on
for k=1:length(data)
    plot(data_in(k,1),data_in(k,2),'k.','MarkerSize',40)
    plot(data_in(k,1),data_in(k,2),'.','Color',data_colors_rgb(k,:),'MarkerSize',35)
end