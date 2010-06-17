%% parameters

% WAV file to analyze
dataFileName = '../data/DDM_speech/WAV/n.wav';

% analysis parameters
numFormants = 2;
numAntiF = 1;
aParams.wType = 'hamming';      % window type
aParams.wLengthMS = 50;         % Length of window (in milliseconds)
aParams.wOverlap = 0.5;         % Factor of overlap of window
aParams.lpcOrder = 12;          % Number of AR coefficients
aParams.zOrder = 4;             % Number of MA coefficients
aParams.peCoeff = 0.9;          % Pre-emphasis factor (0.9, 0.7, etc.)
aParams.fs = 3000;              % sampling rate (in Hz) to resample to

% tracker parameters
cepOrder = 20;
trackBW = 0;
plot_flag = 0; % plot transfer functions
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

% initial state
x0 = [200 980 1000]'; % F Z 
bwStates = [50 50 50]'; % Fbw Zbw

%%
[x_est, x_errVar, x] = runWaveZ(cepOrder, numFormants, numAntiF, ...
    trackBW, dataFileName, algFlag, x0, bwStates, aParams);

%% Super-impose over a spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2(x, x_est, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
axis tight