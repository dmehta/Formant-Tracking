%% runSynth_ARMApq_wrapper.m; DOES NOT TRACK BW
% set parameters for runSynth_ARMApq.m
% 
% INPUT:
%    F:         center frequencies of the resonances (col vector), in Hz
%    Fbw:       corresponding bandwidths of the resonances (col vector), in Hz
%    Z:         center frequencies of the anti-resonances (col vector), in Hz
%    Zbw:       corresponding bandwidths of the anti-resonances (col vector), in Hz
%    dur:       duration of signal, in s
%    pNoiseVar: process noise variance
%    snr_dB:    observation noise, in dB
%    cepOrder:  Number of cepstal coefficients to compute
%    fs:        sampling rate of waveform, in Hz
%    trackBW:   track bandwidths if 1
%    plot_flag: plot figures if 1
%    algFlag:   select 1 to run, 0 not to for [EKF EKS]
%    x0:        initial state of formant trackers [F;Z], in Hz
% 
% OUTPUT:
%    rmse_mean: average RMSE across all tracks

clear 

%% parameters
F = [500]';
Fbw = [100]';
Z = []'; Zbw = []';

dur = .25; % in s
pNoiseVar = 100;
snr_dB = 4;
cepOrder = 15;
fs = 10e3;
trackBW = 0;
plot_flag = 1;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

if trackBW
    x0 = [F; Fbw; Z; Zbw]+100;
else
    x0 = [F; Z]+100;
end

[rmse, x_est] = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
        cepOrder, fs, trackBW, plot_flag, algFlag, x0);
    
figure, hold on
plot(x_est(1:length(F), :)', 'b')
plot(x_est(length(F)+1:end, :)', 'r')
xlabel('Frame')
ylabel('Frequency (Hz)')
title('Estimated trajectories')
rmse