%% synthesize waveforms with VTR trajectories as ground truth, also
% synthesize f0 during voiced frames
% if F3 is greater than 3500, reflect so those values are less than 3500 Hz
% updated 12/16/11 to synthesize f0 correctly

% addpath(genpath('../')); % Paths

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
plotVad = 0;

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
    vtrDbNum = ii;
    %%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load ../data/VTR_Timit/allVTRTracks;
    curData = DATA{vtrDbNum};
    clear DATA

    dataFileName = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.wav');
    dataFileNameOut = strcat('../data/VTRsynthf0 - F3 under 3500/VTRsynth',num2str(vtrDbNum),'.wav');
    
    % Rescale appropriately
    trueStateVTR = 1000*curData.vtrData(:,1:numFormants)';    % Scaling from units
    bwDataVTR    = 500*curData.vtrData(:,5:(4+numFormants))'; % that stored data was in, why not *1000? ******************

    % f0 for voiced frames
    f0FileName   = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.f0');
    data = load(f0FileName);
    f0 = data(:,1)';
    vInd = data(:,2)';
    clear data
    
    % Read in TIMIT audio waveform
    [wav_orig, fs_in] = wavread(dataFileName);
    % Resample input data if sampling rate is not equal to that of input
    if(fs ~= fs_in)
        wav = resample(wav_orig,fs,fs_in,2048);
        %display(['Input Fs = ' int2str(fs_in) ' Hz; resampling to ' num2str(fs) ' Hz'])
    else
        wav = wav_orig;
    end
    
    % Compute window length in samples, now that sampling rate is set
    wLength = floor(aParams.wLengthMS/1000*fs);
    aParams.wLength = wLength + (mod(wLength,2)~=0); % Force even
    
    % simple VAD for speech/silence frames
    %[frameInds] = simpleVAD(wav, aParams, quantThresh, multiBand, [], plotVad);

    % TIMIT labels to help define speech/silence frames
    load([dataFileName(1:end-4), '.mat']) % 0.44 s
    frameInds = timitVAD(wav, data.phones, aParams, plotVad); % 0.3 s

    numFrames = min([size(f0, 2), size(vInd, 2), size(trueStateVTR, 2), size(frameInds, 2)]);

    trueStateVTR = trueStateVTR(1:numFormants,1:numFrames);
    bwDataVTR    = bwDataVTR(1:numFormants,1:numFrames);
    f0           = f0(:,1:numFrames);
    vInd         = vInd(:,1:numFrames);
    frameInds    = frameInds(:,1:numFrames);
    
    F = trueStateVTR;
    
    % 6/20/11 modification to reflect F3 under 3500 Hz
    cutoff = 3500;
    F3 = F(3, 1:end);
    [temp, clipindex] = find(F3 > cutoff);
    F3clip = F3(clipindex);
    F3noclip = cutoff - (F3clip - cutoff);
    F3(clipindex) = F3noclip;
    F(3, 1:end) = F3; % modify formant matrix
    
    Fbw = bwDataVTR;
    f0info = [frameInds; vInd; f0];
    
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
            fprintf([num2str(ii), ' ']) % print database number

            % comment if no variable parameter
            %cepOrder = [cepOrder(1), variableParam(ii)];
            snr_dB = [snr_dB(1), variableParam(kk)];

            [rmse{jj,kk}, x_est{jj,kk}, x_errVar{jj,kk}, x{jj,kk}, trueState{jj,kk}] = runSynth_OLA(F, Fbw, Z, Zbw, N, snr_dB, ...
                    cepOrder, fs, trackBW, plot_flag, algFlag, x0, aParams, sParams, f0info);

            % save synthesized waveform to file
            wavwritedm(x{jj,kk},fs,16,dataFileNameOut)
            %soundsc(x{jj,kk},fs)
            %pause(0.1)
        end
    end
end