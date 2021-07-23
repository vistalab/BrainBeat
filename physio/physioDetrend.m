function tseries_reg = physioDetrend(tseries, dummy_cycles,polyOrder)
% Detrend all the time series elements by removing a 2nd order polynomial
%
% Inputs
%   tseries
%   dummy_cycles
%   polyOrder -  Degree of polynomial we are removing.  Default is 2
%
% Optional key/val
%
% Outputs
%
%
% See also

if notDefined('polyOrder'), polyOrder = 2; end

sliceList =  1:size(tseries,3);
tpoints = (dummy_cycles : nVolumes);
signal_raw_v = single(reshape(tseries(:,:,sliceList,tpoints),[],numel(tpoints)));
signal_reg_v = zeros(size(signal_raw_v));

for voxel = 1:size(signal_raw_v,1)
    signal = signal_raw_v(voxel,:);
    coeff = polyfit(t_vols(tpoints), signal, polyOrder);
    signal_fit = polyval(coeff, t_vols(tpoints));
    
    % We work on the residual
    signal_reg_v(voxel,:) = signal - signal_fit;
end

tseries_reg = reshape(signal_reg_v, size(tseries(:,:,sliceList,tpoints)));

end