% synthesize waveforms with VTR trajectories as ground truth, stochastic input at each frame
% addpath(genpath('../')); % Paths
% if F3 is greater than 3500, reflect so those values are less than 3500 Hz

%%
clear all

% synthesis parameters
% dur = 2; % in s
fs = 16000; % in Hz
sParams.wType = 'hanning';      % window type
sParams.wLengthMS  = 20;       % Length of window (in milliseconds)
sParams.wOverlap = 0.5;         % Factor of overlap of window
snr_dB(1) = 20;
cepOrder(1) = 20; %max(aParams.lpcOrder, aParams.zOrder);
randn('state', 5)

numFormants = 4;
Z = []; Zbw = [];

% analysis parameters
aParams.wType = 'hamming';          % window type
aParams.wLengthMS = sParams.wLengthMS;  % Length of window (in milliseconds)
aParams.wOverlap = sParams.wOverlap;    % Factor of overlap of window
aParams.lpcOrder = 12;      % Number of AR coefficients
aParams.zOrder = 0;        % Number of MA coefficients
aParams.peCoeff = 0.7;                    % Pre-emphasis factor (0.9, 0.7, etc.)
aParams.fs = fs;                        % sampling rate (in Hz) to resample to

% tracker parameters
snr_dB(2) = 20;
cepOrder(2) = 20; %max(aParams.lpcOrder, aParams.zOrder);
trackBW = 0;
plot_flag = 0; % plot transfer functions
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
offset = 100; % set initial state offset, in Hz
variableParam = 20;

% initial state
initFormant = 500 + 1000*(0:(numFormants - 1))';
initBW      = 80 + 40*(0:(numFormants - 1))';
if(trackBW)
    x0 = [initFormant', initBW']';
else
    x0 = initFormant;
end

% Monte Carlo analysis parameters
numTrials = 1;

for ii = 1:516
    %%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vtrDbNum = ii;
    load ../data/VTR_Timit/allVTRTracks;
    curData = DATA{vtrDbNum};

    % Rescale appropriately
    trueStateVTR = 1000*curData.vtrData(:,1:numFormants)';    % Scaling from units
    bwDataVTR    = 500*curData.vtrData(:,5:(4+numFormants))'; % that stored data was in, why not *1000? bec Dan said Abeer said so******************

    F = trueStateVTR;
    Fbw = bwDataVTR;
    numFrames = size(F,2);

    % 6/20/11 modification to reflect F3 under 3500 Hz
    cutoff = 3500;
    F3 = F(3, 1:end);
    [temp, clipindex] = find(F3 > cutoff);
    F3clip = F3(clipindex);
    F3noclip = cutoff - (F3clip - cutoff);
    F3(clipindex) = F3noclip;
    F(3, 1:end) = F3; % modify formant matrix
    
    wLength = floor(sParams.wLengthMS/1000*fs);
    wLength = wLength + (mod(wLength,2)~=0); % Force even
    win = feval(aParams.wType,wLength);
    N = (numFrames+1)*wLength*(1-sParams.wOverlap);
    
    %% synthesize and track
    rmse = cell(numTrials, 1);
    x_est = cell(numTrials, 1);
    x_errVar = cell(numTrials, 1);
    x = cell(numTrials, 1);
    trueState = cell(numTrials, 1);

    % loop through
    for kk = 1:length(variableParam)
        for jj = 1:numTrials
            %disp(['Processing Trial #', num2str(jj), '...'])
            fprintf('.')

            % comment if no variable parameter
            %cepOrder = [cepOrder(1), variableParam(ii)];
            snr_dB = [snr_dB(1), variableParam(kk)];

            [rmse{jj,kk}, x_est{jj,kk}, x_errVar{jj,kk}, x{jj,kk}, trueState{jj,kk}] = runSynth_OLA(F, Fbw, Z, Zbw, N, snr_dB, ...
                    cepOrder, fs, trackBW, plot_flag, algFlag, x0, aParams, sParams);

            % save synthesized waveform to file
            synthFileName   = strcat('../data/VTRsynth - F3 under 3500/VTRsynth',num2str(vtrDbNum),'.wav');
            wavwritedm(x{jj,kk},fs,16,synthFileName)
        end
    end
end