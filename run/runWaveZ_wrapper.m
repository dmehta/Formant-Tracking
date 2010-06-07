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
cepOrder = 25;
numFormants = 3;
numAntiF = 2;
trackBW = 0;
% dataFileName = '../data/synthData/mlm.tea.10k.wav';
% dataFileName = '../data/DDM_speech/WAV/ah.wav';
dataFileName = '../data/DDM_speech/WAV/n.wav';
% dataFileName = '../data/DDM_speech/WAV/m.wav';
% dataFileName = '../data/DDM_speech/WAV/ana.wav';
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

%% initial state
% initFormant = 200 + 1000*(0:(numFormants - 1))';
% initBW = 80 + 40*(0:(numFormants - 1))';
% 
% initAntiF = 1200 + 1000*(0:(numAntiF - 1))';
% initAntiFBW = 80 + 40*(0:(numAntiF - 1))';

initFormant = [200 2000 3500]';
initBW = [100 100 100]';

initAntiF = [600 2500]';
initAntiFBW = [100 100]';

if numFormants && numAntiF
    if trackBW
        x0 = [initFormant; initBW; initAntiF; initAntiFBW];
    else
        x0 = [initFormant; initAntiF];
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
        x0 = [initAntiF; initAntiFBW];
    else
        x0 = [initAntiF];
    end
end

%%
tic
x_est = runWaveZ(cepOrder, numFormants, numAntiF, ...
    trackBW, dataFileName, algFlag, x0);
toc

%% plot all tracks
figure, hold on
for jj = 1:numTrials
    plot(x_est(1:numFormants, :)', 'b')
    plot(x_est(numFormants+1:end, :)', 'r')
end
xlabel('Frame')
ylabel('Frequency (Hz)')
title('Estimated EKS trajectories (Formant-blue, Antiformant-red)')
    