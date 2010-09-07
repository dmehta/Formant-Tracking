%% parameters

% WAV file to analyze
% dataFileName = '../data/DDM_speech/n.wav';
dataFileName = '../data/DDM_speech/ana.wav';
% dataFileName = '../data/synthData/ln.roy.10k.wav';
% dataFileName = '../data/synthData/tfq.roy.10k.wav';

% analysis parameters
numFormants = 3;
numAntiF = 0;
aParams.wType = 'hamming';      % window type
aParams.wLengthMS = 20;         % Length of window (in milliseconds)
aParams.wOverlap = 0.5;         % Factor of overlap of window
aParams.lpcOrder = 16;          % Number of AR coefficients
aParams.zOrder = 0;             % Number of MA coefficients
aParams.peCoeff = 0;          % Pre-emphasis factor (0.9, 0.7, etc.)
aParams.fs = 8000;              % sampling rate (in Hz) to resample to

% tracker parameters
cepOrder = 16; %max(aParams.lpcOrder, aParams.zOrder);
trackBW = 1;
plot_flag = 0; % plot transfer functions
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

% initial state
% x0 = [500 1000 1500 1500 3541]'; % F Z 
% bwStates = [50 50 50 50 50]'; % Fbw Zbw

% x0 = [500 1000 1500 1100]'; % F Z 
% bwStates = [50 50 50 50]'; % Fbw Zbw

x0 = [500 1500 2500]'; % F Z 
bwStates = [50 50 50]'; % Fbw Zbw

%%
[x_est, x_errVar, x] = runWaveZ(cepOrder, numFormants, numAntiF, ...
    trackBW, dataFileName, algFlag, x0, bwStates, aParams);

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
set(gca, 'PlotBoxAspectRatio', [3 1 1])

%% Praat tracks on spectrogram
praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
params = getFormantParamsDefault;
 params.winLen = aParams.wLengthMS/1000;
 params.timestep = params.winLen*aParams.wOverlap;
 params.numTracks = numFormants;
[f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(x, aParams.fs, praat_dir, params);
plotSpecTracksWS(x, f', bw', aParams);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
format_plot
set(gca, 'PlotBoxAspectRatio', [3 1 1])
