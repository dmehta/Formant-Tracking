%% runSynth_wrapper.m
% set parameters for runSynth.m
% 
% INPUT
%
% testMethod    : Synth, VTR, PRAAT, WS
% snr_dB        : How much observation noise to add
% cepOrder     : How many cepstal coefficients to include in the observations
% fs            : Sampling rate at which the observations are made (fake here)
% numFormants   : Number of formants that we should track
% trackBW       : Flag whether to track (1) or not track (0) bandwidths
% numParticles  : ?
% doPlots       : Flag for whether to plot (1) or not to plot (0) estimates
%                       and ground truth tracks
% varargin      : Dependent on testMethod. If
%                       Synth: number of observations, process noise variance
%                       VTR: data filename, sample index
%                       PRAAT, WS: data filename

%% scenario 1
testMethod = 'Synth';
snr_dB = 100;
cepOrder = 15;
fs = 16e3;
numFormants = 4;
trackBW = 1;
numParticles = [];
N = 30;
pNoiseVar = 100;
doPlots = 1;

rmse = runSynth(testMethod, snr_dB, cepOrder, fs, numFormants, trackBW, numParticles, doPlots, N, pNoiseVar);

%%
