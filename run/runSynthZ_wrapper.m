%% runSynthZ_wrapper.m
% set parameters for runSynthZ.m
% 
% INPUT
%
% nP  : Number of pole pairs to track
% nZ  : Number of zero pairs to track
% N    : Number of observations to generate
% pNoiseVar : Process noise variance
% snr_dB    : How much observation noise to add
% cepOrder : How many cepstal coefficients to include in the observations
% fs        : Sampling rate at which the observations are made (fake here)

%% scenario 1
nP = 4;
nZ = 4;
N = 30;
pNoiseVar = 1;
snr_dB = 100;
cepOrder = 15;
fs = 16e3;
%trackBW = 0; % don't know how to track BW yet
%doPlots = 1;

rmse = runSynthZ(nP, nZ, N, pNoiseVar, snr_dB, cepOrder, fs);

%%
