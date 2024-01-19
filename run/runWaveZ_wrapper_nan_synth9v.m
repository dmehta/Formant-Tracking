% addpath(genpath('../'));

%% parameters

clear all

% WAV file to analyze
% dataFileName = '../data/DDM_speech/n.wav';
% dataFileName = '../data/DDM_speech/ana.wav';
% dataFileName = '../data/synthData/ln.roy.10k.wav';
% dataFileName = '../data/synthData/tfq.roy.10k.wav';
% dataFileName = '../data/synthData/nan_synth2.wav';
% dataFileName = '../data/synthData/nan_synth4.wav';
% dataFileName = '../data/synthData/nan_synth5.wav'; % stochastic input
% dataFileName = '../data/synthData/nan_synth6.wav'; % periodic input
% dataFileName = '../data/synthData/nan_synth7.wav'; % stochastic input
% dataFileName = '../data/synthData/nan_synth8.wav'; % stochastic input,
% 800 Hz zero
dataFileName = '../data/synthData/nan_synth9.wav'; % stochastic input, different /n/2
% dataFileName = '../data/synthData/nan_synth9v.wav'; % periodic input, different /n/2
% dataFileName = '../data/synthData/nan_synth10.wav'; % stochastic input, n2 with 20ms seg

% analysis parameters
aParams.wType = 'hamming';      % window type
aParams.wLengthMS = 100;         % Length of window (in milliseconds)
aParams.wOverlap = 0.5;         % Factor of overlap of window
aParams.lpcOrder = 4;          % Number of AR coefficients
aParams.zOrder = 2;             % Number of MA coefficients
aParams.peCoeff = 0.9;          % Pre-emphasis factor (0.9, 0.7, etc.)
aParams.fs = 10000;              % sampling rate (in Hz) to resample to

% tracker parameters
numFormants = 2;
numAntiF = 1;
cepOrder = 20; %max(aParams.lpcOrder, aParams.zOrder);
trackBW = 1;
plot_flag = 0; % plot transfer functions
algFlag = [1 0]; % Select 1 to run, 0 not to; [EKF EKS] NOT READY FOR EKS YET

formantInds = ones(2*numFormants+2*numAntiF, 75); % track everything, f/bw
formantInds(5:6, 20:56) = 0; % coasting antiformant f and bw

% formantInds = ones(numFormants+numAntiF, 75); % track everything, f
% formantInds(3, 20:56) = 0; % coasting antiformant, not tracking bw

% initial state
% x0 = [500 1000 1500 1500 3541]'; % F Z 
% bwStates = [50 50 50 50 50]'; % Fbw Zbw

% x0 = [500 1000 1500 1100]'; % F Z 
% bwStates = [50 50 50 50]'; % Fbw Zbw

% x0 = [500 1500 2500]'; % F Z 
% bwStates = [50 50 50]'; % Fbw Zbw

x0 = [500 1500 1000]'; % F Z 
bwStates = [32 100 52]'; % Fbw Zbw

%%
[x_est, x_errVar, x] = runWaveZ_nan_synth9v(cepOrder, numFormants, numAntiF, ...
    trackBW, dataFileName, algFlag, x0, bwStates, aParams, formantInds);

%% Super-impose over a spectrogram
load trueState9v
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

% plotSpecTracks2BW(x, x_est, aParams, numAntiF, trackBW); % with bandfwidths
% plotSpecTracks2_estVar(x, x_est, x_errVar, aParams, numAntiF, trackBW); % with uncertainty
plotSpecTracks2_estVar_Truth(x, x_est, x_errVar, trueState, aParams, numAntiF, trackBW); % with truth and uncertainty

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
plotSpecTracksPraat(x, f', bw', aParams, formantTimes);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
format_plot
set(gca, 'PlotBoxAspectRatio', [3 1 1])

%% tracks with estimated variances and truth
% trueState = 1000*ones(size(x_est));
load trueState9v

titleCell(1,2) = {'EKS'}; % hard coded for now
titleCell(2,2) = {'b:'};
titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};      % Color for true state plot
nP = numFormants;

plotStateTracksFZ_EstVar_Truth(trueState,x_est,x_errVar,titleCell,nP,trackBW)

% Compute and Display RMSE
rmse = zeros(size(x_est, 1),1);
for j = 1:size(x_est, 1)
    meas = x_est(j,:);
    truest = trueState(j,:);
    
    meas(find(isnan(meas))) = [];
    truest(find(isnan(truest))) = [];
    
    rmse(j) = norm(meas-truest)/sqrt(size(x_est,2));
end

rmse
mean(rmse)