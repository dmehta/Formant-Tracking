%% runSynth_ARMApq_wrapper.m
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
%    plot_flag: plot figures if 1
%    algFlag:   select 1 to run, 0 not to for [EKF EKS]
%    x0:        initial state of formant trackers [F;Z], in Hz
% 
% OUTPUT:
%    rmse_mean: average RMSE across all tracks

%% parameters
F = [500 1000]';
Fbw = [100 200]';
Z = []; Zbw = [];

% F = [1500]';
% Fbw = [100]';
% Z = [1000 2000]'; Zbw = [10 50]';

dur = .5; % in s
pNoiseVar = 10;
snr_dB = 25;
cepOrder = 15;
fs = 16e3;
plot_flag = 1;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
x0 = [F; Z]+0;

[rmse_EKF, x_estEKF] = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
        cepOrder, fs, plot_flag, algFlag, x0);
% [rmse_EKF, x_estEKF, rmse_EKS, x_estEKS] = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, ...
%         cepOrder, fs, plot_flag, algFlag, x0);
    
figure, hold on
plot(x_estEKF(1:length(F), :)', 'b')
plot(x_estEKF(length(F)+1:end, :)', 'r')
xlabel('Frame')
ylabel('Frequency (Hz)')
title('Estimated EKF trajectories')

% figure, hold on
% plot(x_estEKS(1:length(F), :)', 'b')
% plot(x_estEKS(length(F)+1:end, :)', 'r')
% xlabel('Frame')
% ylabel('Frequency (Hz)')
% title('Estimated EKS trajectories')