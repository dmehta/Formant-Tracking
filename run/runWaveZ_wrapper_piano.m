% addpath(genpath('../'));

%% parameters

clear all

% WAV file to analyze
dataFileName = '../data/DDM_speech/piano.wav';

% analysis parameters
aParams.wType = 'hamming';      % window type
aParams.wLengthMS = 20;         % Length of window (in milliseconds)
aParams.wOverlap = 0.5;         % Factor of overlap of window
aParams.lpcOrder = 16;          % Number of AR coefficients
aParams.zOrder = 4;             % Number of MA coefficients
aParams.peCoeff = 0.7;          % Pre-emphasis factor (0.9, 0.7, etc.)
aParams.fs = 8000;              % sampling rate (in Hz) to resample to

% tracker parameters
numFormants = 3;
numAntiF = 2;
cepOrder = 16; %max(aParams.lpcOrder, aParams.zOrder);
trackBW = 1;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

formantInds = []; % don't use speech detection
% formantInds = 1; % use speech detection

% formantInds = ones(2*numFormants+2*numAntiF, 787); % track everything, f/bw

% coasting antiformant f and bw, uncomment 2 lines:
% formantInds(2*numFormants+1:end, [1:451 582:end]) = 0; % coasting antiformant f and bw
% algFlag = [1 0]; % Select 1 to run, 0 not to; [EKF EKS]; EKS not coded for changing state vector; only EKF 

% initial state
% x0 = [500 1500 2500]'; % F Z 
% bwStates = [50 50 50]'; % Fbw Zbw

x0 = [500 1500 2500 1000 2000]'; % F Z 
bwStates = [80 120 160 80 80]'; % Fbw Zbw

%%
% to keep state vector length the same
[x_est, x_errVar, x] = runWaveZ(cepOrder, numFormants, numAntiF, ...
    trackBW, dataFileName, algFlag, x0, bwStates, aParams, formantInds);

% to change state vector length
% [x_est, x_errVar, x] = runWaveZ_piano(cepOrder, numFormants, numAntiF, ...
%     trackBW, dataFileName, algFlag, x0, bwStates, aParams, formantInds);

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