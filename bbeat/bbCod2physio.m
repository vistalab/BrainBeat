function [out_r_map]=bbCod2physio(ni_odd,ni_even,slices)
%
% Function calculates coefficient of determination between PPG triggered
% signal for even and odd heartbeats.
%
% [out_r_map] = bbCod2physio(ni)
%
% required input:
% ni: a nifti structure loaded by niftiRead
%
% optional input:
% slices: slices for which to calculate response function, defaults = do
%           all slices
% 
% output:
% out_r_map: voxels X voxels X slices maps with 'correlation' to hearbeat
%
% Wandell Copyright Vistasoft Team, 2013
% Written by Dora Hermes 2014

if ~exist('slices','var') % do whole brain
    slices = [1:size(ni_odd.data,3)];
end

% create output matrix:
out_r_map = zeros(size(ni_odd.data,1),size(ni_odd.data,2),length(slices));

for s = 1:length(slices)
    disp(['slice ' int2str(s) ' of ' int2str(length(slices))])
    sli = slices(s);

    %%%%% calculate correlation map
    odd_resp = squeeze(ni_odd.data(:,:,sli,:));
    odd_resp = reshape(odd_resp,[size(odd_resp,1) * size(odd_resp,2),size(odd_resp,3)]);
    
    even_resp = squeeze(ni_even.data(:,:,sli,:));
    even_resp = reshape(even_resp,[size(even_resp,1) * size(even_resp,2),size(even_resp,3)]);
    
    r = calccod(odd_resp',even_resp',1,0,0)./100; 
    % make sure values go between 0 and 1, smaller than 0 is just worse...
    r(r<0) = 0;
    r(r>1) = 1;
    r = reshape(r,[size(ni_odd.data,1),size(ni_odd.data,2)]);
    
    out_r_map(:,:,s) = r;
    
    clear r odd_resp even_resp
    
end
