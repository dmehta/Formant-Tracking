% addpath(genpath('../'));

%% parameters

% clear all

% WAV file to analyze
% dataFileName = 'temp5.wav';

dataFileName = '../data/DDM_speech/n.wav';
% dataFileName = '../data/DDM_speech/ana.wav';
% dataFileName = '../data/DDM_speech/m.wav';
% dataFileName = '../data/DDM_speech/ng.wav';
% dataFileName = '../data/DDM_speech/a.wav';
% dataFileName = '../data/DDM_speech/a_wh.wav';
% dataFileName = '../data/DDM_speech/piano.wav';

% dataFileName = '../data/synthData/ln.roy.10k.wav';
% dataFileName = '../data/synthData/tfq.roy.10k.wav';

% dataFileName = '../data/synthData/nan_synth2.wav';
% dataFileName = '../data/synthData/nan_synth4.wav';
% dataFileName = '../data/synthData/nan_synth5.wav'; % stochastic input
% dataFileName = '../data/synthData/nan_synth6.wav'; % periodic input
% dataFileName = '../data/synthData/nan_synth7.wav'; % stochastic input
% dataFileName = '../data/synthData/nan_synth8.wav'; % stochastic input,
% 800 Hz zero
% dataFileName = '../data/synthData/nan_synth9.wav'; % stochastic input,
% different /n/2
% dataFileName = '../data/synthData/nan_synth10.wav'; % stochastic input, n2 with 20ms seg

vtrDbNum = 1;
dataFileName = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.wav');
% dataFileName = strcat('../data/VTRsynth/VTRsynth',num2str(vtrDbNum),'.wav');

% analysis parameters
aParams.wType = 'hamming';      % window type
aParams.wLengthMS = 20;         % Length of window (in milliseconds)
aParams.wOverlap = 0.5;         % Factor of overlap of window
aParams.lpcOrder = 12;          % Number of AR coefficients
aParams.zOrder = 4;             % Number of MA coefficients
aParams.peCoeff = 0.7;          % Pre-emphasis factor (0.9, 0.7, etc.)
aParams.fs = 7000;              % sampling rate (in Hz) to resample to

% tracker parameters
numFormants = 3;
numAntiF = 1;
cepOrder = 15; %max(aParams.lpcOrder, aParams.zOrder);
trackBW = 1;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

% formantInds = ones(2*numFormants+2*numAntiF, 75); % track everything, f/bw
% formantInds(5:6, 20:56) = 0; % coasting antiformant f and bw

% formantInds = ones(numFormants+numAntiF, 75); % track everything, f
% formantInds(3, 20:56) = 0; % coasting antiformant, not tracking bw

formantInds = []; % don't use speech detection
% formantInds = 1; % use speech detection

% initial state
% x0 = [500 1000 1500 1500 3541]'; % F Z 
% bwStates = [50 50 50 50 50]'; % Fbw Zbw

% x0 = [500 1000 1500 1100]'; % F Z 
% bwStates = [50 50 50 50]'; % Fbw Zbw

% x0 = [500 1500 2500]'; % F Z 
% bwStates = [50 50 50]'; % Fbw Zbw

% x0 = [500 1500 2500 1000 2000]'; % F Z 
% bwStates = [80 120 160 80 80]'; % Fbw Zbw

x0 = [500 1500 2500 1000]'; % F Z 
bwStates = [80 120 160 80]'; % Fbw Zbw

% x0 = [500 2000 2510 3400]'; % F Z 
% bwStates = [80 80 80 80]'; % Fbw Zbw

%%
[x_est, x_errVar, x] = runWaveZ(cepOrder, numFormants, numAntiF, ...
    trackBW, dataFileName, algFlag, x0, bwStates, aParams, formantInds);

%% Super-impose over a spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2BW(x, x_est, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
format_plot
set(gca, 'PlotBoxAspectRatio', [3 1 1])

mean(x_est,2)

%% tracks with estimated variance
plotStateTracksFZ_EstVar(x_est,x_errVar,numFormants,trackBW)

%% Wavesurfer tracks on spectrogram
[f,bw] = wavesurferFormantRead([dataFileName(1:end-3) 'frm'],numFormants);
plotSpecTracksWS(x, f, bw, aParams)
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
format_plot
% set(gca, 'PlotBoxAspectRatio', [3 1 1])

%% Praat tracks on spectrogram
praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
params = getFormantParamsDefault;
 params.winLen = aParams.wLengthMS/1000;
 params.timestep = params.winLen*aParams.wOverlap;
 params.numTracks = numFormants;
[f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(x, aParams.fs, praat_dir, params);
plotSpecTracksPraat(x, f', bw', aParams, formantTimes);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
format_plot
% set(gca, 'PlotBoxAspectRatio', [3 1 1])

%% EKS with uncertainty on spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2_estVar(x, x_est, x_errVar, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
set(gca, 'PlotBoxAspectRatio', [3 1 1])