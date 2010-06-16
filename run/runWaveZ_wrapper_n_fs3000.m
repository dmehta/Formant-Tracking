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
cepOrder = 10;
numFormants = 2;
numAntiF = 1;
trackBW = 0;
dataFileName = '../data/DDM_speech/WAV/n.wav';
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]

%% initial state
initFormant = [213.57115      981.71228]';
initBW = [40.56219      7.451091]';

initAntiF = [1000]';
initAntiFBW = [50]';

%%
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
[x_est, aParams] = runWaveZ(cepOrder, numFormants, numAntiF, ...
    trackBW, dataFileName, algFlag, x0);

%% Super-impose over a spectrogram
[x, fs] = wavread(dataFileName);
x = resample(x,3000,fs,2048);
plotSpecTracks2(x, x_est, aParams, numAntiF, trackBW);
axis tight
% ylim([0 3000])