function sig_conv = physioResponse(hbeats, tseries, TR, nvols, slice_timing, coor, dummy_cycles)

if ~exist('dummy_cycles', 'var')
    dummy_cycles = 3;
end

t_vols = (0 : (nvols - 1)) * TR;

% regress the fmri time series 
slices = 1:size(tseries,3);
tpoints = dummy_cycles : size(tseries,4);
signal_raw_v = single(reshape(tseries(:,:,slices,tpoints),[],numel(tpoints)));
signal_reg_v = zeros(size(signal_raw_v));
for voxel = 1:size(signal_raw_v,1)
    signal = signal_raw_v(voxel,:);
    coeff = polyfit(t_vols(tpoints), signal, 2);
    signal_fit = polyval(coeff, t_vols(tpoints));
    signal_reg_v(voxel,:) = signal - signal_fit;
end

tseries_reg = reshape(signal_reg_v, size(tseries(:,:,slices,tpoints)));

% find the fmri volume closest to each ppg trigger
tt = t_vols(dummy_cycles:end) + slice_timing(sl);    % correct for sliceTiming for each slice
hbeats = hbeats(hbeats > TR * (dummy_cycles-1));
for trig = 1:numel(hbeats)
    tseries_vol_idx = find(tt > hbeats(trig), 1, 'first');
    tseries_vol_t = tt(tseries_vol_idx);
    if ~isempty(tseries_vol_idx)
        delta_t(trig) = tseries_vol_t - hbeats(trig);
%         if delta_t(trig) > 2000
%             keyboard;
%         end
        sig_hbeats(trig) = tseries_reg(x, y , sl, tseries_vol_idx);
    end
end

%figure;plot(delta_t, fmri_sig, 'x')
    
[delta_t_sort, idx] = sort(delta_t);
sig_hbeats_sort = sig_hbeats(idx);

%figure;plot(new_t, new_sig, 'x')
new_window = 10:10:2000;
sig_interp = interp1(delta_t_sort+rand(size(delta_t_sort)), sig_hbeats_sort, new_window);
%figure;plot((10:10:2000), interp_sig,'x')
%ker = fspecial('gaussian', [1,8]);
ker = ones(1, 16);
sig_conv = conv(sig_interp, ker, 'same');


end

