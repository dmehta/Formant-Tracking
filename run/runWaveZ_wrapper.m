% set parameters for runWaveZ.m

% INPUT:
%    cep_order:     How many cepstal coefficients to include in the observations
%    fs_in:         Sampling rate, in Hz
%    numFormants:   Number of formants (pole pairs) that we should track
%    numAntiF:      Number of anti-formants (zero pairs) that we should track
%    trackBW:       Flag whether to track (1) or not track (0) bandwidths;
%                       right now, must be set to 1 because we do not know how to set
%                       bandwidths of anti-formants
%    dataFileName:  Audio (.wav) filename
%    wsFileName:    Currently compared against Wavesurfer. This must contain
%                       at least as many formants as we wish to track and
%                       preferably the exact number to avoid mismatch in
%                       resampling issues that might come up.
%    algFlag:       Select 1 to run, 0 not to for [EKF EKS]

clear 

%% parameters
cepOrder = 15;
fs_in = 2e3;
numFormants = 1;
numAntiF = 1;
trackBW = 0;
% dataFileName = '../data/synthData/mlm.tea.10k.wav';
% dataFileName = '../data/DDM_speech/WAV/ah.wav';
dataFileName = '../data/DDM_speech/WAV/n.wav';
wsFileName = '../data/synthData/mlm.tea.10k.FRM';
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

%% initial state
initFormant = 500 + 1000*(0:(numFormants - 1))';
initBW = 80 + 40*(0:(numFormants - 1))';

initAntiF = 700;
initAntiFBW = 80;

if numFormants && numAntiF
    if trackBW
        x0 = [initFormant; initBW; initAntiF; initAntiFBW];
    else
        x0 = [];
    end
end

if numFormants && ~numAntiF
    if trackBW
        x0 = [initFormant; initBW];
    else
        x0 = initFormant;
    end
end

if ~numFormants && numAntiF
    if trackBW
        x0 = [];
    else
        x0 = [];
    end
end

%%
numTrials = 1;
rmse = zeros(numTrials, 1);
x_est = cell(numTrials, 1);

for jj = 1:numTrials
    tic
    [rmse(jj), x_est{jj}] = runWaveZ(cepOrder, fs_in, numFormants, numAntiF, ...
        trackBW, dataFileName, wsFileName, algFlag, x0);
    toc
end

%% plot all tracks
figure, hold on
for jj = 1:numTrials
    plot(x_est{jj}(1:numFormants, :)', 'b')
    plot(x_est{jj}(numFormants+1:end, :)', 'r')
end
xlabel('Frame')
ylabel('Frequency (Hz)')
title('Estimated EKS trajectories (Formant-blue, Antiformant-red)')
    