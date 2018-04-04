s_nr = [2 3 4];
scan_nr  = [1 1 2];

all_pcs = zeros(length(s_nr),2,128);
for kk = 1:length(s_nr)
    load(['./local/s-' int2str(s_nr(kk)) '_scan-' int2str(scan_nr(kk)) 'pc12'],'y1','y2','t_hr')
    all_pcs(kk,1,:) = y1;
    all_pcs(kk,2,:) = y2;
end


%% 

all_pcs = reshape(all_pcs,size(all_pcs,1)*size(all_pcs,2),size(all_pcs,3))';

%%

[u,s,v] = svd(all_pcs);
% s = diag(s);

figure
temp = u*s;
% plot(u(:,1:3))
plot(temp(:,1:3))

% TODO: Save first two, run through all voxels, all subjects...