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
snr_dB = 105;
cepOrder = 25;
fs = 16e3;
numFormants = 2;
trackBW = 1;
numParticles = [];
doPlots = 1;
N = 30;
pNoiseVar = 1000;

%%
rmse = runSynth(testMethod, snr_dB, cepOrder, fs, numFormants, trackBW, numParticles, doPlots, N, pNoiseVar);

%% Dan's Fig. 3.5
runSynth('WS', 15, 15, 10000, 3, 1, [], 1, '../data/synthData/ln.roy.10k.FRM');

%%
runSynth('Synth',105, 15, 16000, 2, 1, [], 1, 30, 1000);