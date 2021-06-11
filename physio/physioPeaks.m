function [onsets, band_sig] = physioPeaks(data, dt)
% Calculate the peaks for either the PPG or RESP data
% See also physioGet.m
%
% Input: 
%   data        physio waveform, ppg or resp
%   dt          time resolution of the waveform
%
% Output:
%   onsets      detected triggers in milliseconds
%   band_sig    bandpass filtered waveform
%

srate = 1000/dt;
% find the peaks in the autocorrelation function:
[peaks_ac,peaks_ac_i] = findpeaks(autocorr(data,'NumLags',500));
% the first maximum peak (zero is not included) is the first
% autocorrelation time
[~,max_peaks_ac_i] = max(peaks_ac);
% set the interval at 70% of the ppg/resp rate:
interval = .7 * peaks_ac_i(max_peaks_ac_i)/srate;

% low-pass filter 
% not necessary, keep in for lower quality data?
band = 5;
Rp   = 3; Rs = 60; % third order Butterworth
high_p =  band(1)*2/srate;
delta = 0.001*2/srate;
high_s = min(1-delta,high_p+0.1);

[n_band,wn_band] = buttord(high_p,high_s,Rp,Rs);
[bf_b,bf_a] = butter(n_band,wn_band,'low');
band_sig    = filtfilt(bf_b,bf_a,data);

% detect peak sample positions (in milliseconds)
[~,onsets] = findpeaks(band_sig,'minpeakdistance',round(interval*srate));
onsets = onsets/srate * 1000;

end
