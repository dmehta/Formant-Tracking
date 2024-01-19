% karma_demo.m
%
% Author: Daryush Mehta
% Created: 07/29/2011
% Last modified: 08/01/2011
% 
% Wrapper function around karma() to illustrate usage of functions related to Kalman-based
% autoregressive moving average modeling and inference for formant and
% antiformant tracking.
% 
% Uncomment one of the numbered cell blocks in the code to see how to call
% karma(). Modify as desired.
% 
% References:
% [1]   D. D. Mehta, D. Rudoy, and P. J. Wolfe, "Kalman-based autoregressive
%       moving average modeling and inference for formant and antiformant
%       tracking," In review, 2011.
% [2]   D. Rudoy, D. N. Spendley, and P. J. Wolfe, "Conditionally linear
%       Gaussian models for estimating vocal tract resonances," Proceedings of
%       Interspeech, Antwerp, Belgium, 2007.
% [3]   D. Rudoy, "Nonstationary time series modeling with application to speech
%       signal processing," Doctor of Philosophy thesis, School of Engineering
%       and Applied Sciences, Harvard University, Cambridge, MA, 2010.
%       Chapter 3.

function karma_demo()

% path check
if ~exist('plotSpecTracks2_estVar.m', 'file')
    addpath(genpath('../')); % executes if directories not on path
end

%% 1. run with default parameters
[x_est, x_errVar, x, params] = karma('\..\data\TIMIT\Timit1.wav');

%% 2. Alternatively, specify select modifications to the default parameter set:
% aParams.lpcOrder = 12;
% cepOrder = 15;
% cepType = 1;
% 
% [x_est, x_errVar, x, params] = karma('\..\data\TIMIT\Timit1.wav', [], [], aParams, cepOrder, cepType);

%% 3. Or, specify all of the parameters:
% numFormants = 3;
% numAntiF = 0;
% aParams.peCoeff = 0.7;
% aParams.wType = 'hamming';
% aParams.wLengthMS = 20;
% aParams.wOverlap = 0.5;
% aParams.lpcOrder = 12;
% aParams.zOrder = 0;
% aParams.fs = 7000;
% cepOrder = 15;
% cepType = 1;
% algFlag = 2;
% 
% x0 = 500 + 1000*(0:(numFormants - 1))';
% x02 = 1000 + 1000*(0:(numAntiF - 1))';
% x0 = [x0; x02];
% 
% bwStates = 80 + 40*(0:(numFormants - 1))';
% bwStates2 = 80 + 40*(0:(numAntiF - 1))';
% bwStates = [bwStates; bwStates2];
% 
% formantInds = 1;
% trackBW = 1;
% 
% [x_est, x_errVar, x, params] = karma('\..\data\TIMIT\Timit1.wav', numFormants, numAntiF, ...
%     aParams, cepOrder, cepType, algFlag, x0, bwStates, formantInds, trackBW);

%% 4. Here is an antiformant tracking demo with synthesized waveform 
%% (cf. Fig. 4 in [1])
% numFormants = 2;
% numAntiF = 1;
% aParams.lpcOrder = 16;
% aParams.zOrder = 2;
% aParams.fs = 10000;
% cepOrder = 16;
% disp('A few minutes needed to process due to ARMA model fitting...')
% [x_est, x_errVar, x, params] = karma('\..\data\synthData\nan_synth9v.wav', numFormants, numAntiF, aParams, cepOrder);

%% plot routines
% a. Frequency and bandwidth tracks with estimated variance
% +/- 1 standard deviation as thickness of tracks
plotStateTracksFZ_EstVar(x_est,x_errVar,params.numFormants,params.trackBW)

% b. Superimpose tracks on spectrogram with bandwidths
% +/- 3 dB bandwidth as thickness of tracks
plotSpecTracks2BW(x, x_est, params.aParams, params.numAntiF, params.trackBW);
set(gca, 'PlotBoxAspectRatio', [5 1 1])

% c. Superimpose tracks on spectrogram with uncertainty
% +/- 1 standard deviation as thickness of tracks
plotSpecTracks2_estVar(x, x_est, x_errVar, params.aParams, params.numAntiF, params.trackBW);
set(gca, 'PlotBoxAspectRatio', [5 1 1])