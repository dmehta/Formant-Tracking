function [] = runWave(cepOrder, fs_in, numFormants, trackBW, useVAD, useCorr, dataFileName, wsFileName)

% Test tracking algorithms using raw audio waveform. Requires a
% post-processed file from WaveSurfer in order to input bandwidths,
% estimate inter-formant correlation (if desired) and validate the results.
% Nothing here depends on the VTR database, if you wish to use the VTR
% database for testing use testVTR.m
%
% Author: Daniel Rudoy
% Created:  12/13/2007
% Modified: 04/09/2008, 01/07/2007, 02/10/2010

% INPUT
%
% cep_order  :   How many cepstal coefficients to include in the observations
% fs_in      :   Sampling rate at which the observations are made (fake here)
% numFormants:   Number of formants that we should track
% useVAD     :   Use voice activity detection? ('1' - yes, '0' - no)
% useCorr    :   Use correlation? ('1' - yes, '0' - no)
% dataFileName: Audio (.wav) source file
% wsFileName  : Tracks of formants and bandwidths produced by Wavesurfer which
%               is necessary for bootsrapping the inter-formant correlation procedure
%               This file must contain at least as many formants as we wish to track and
%               preferrably the exact number to avoid mismatch in resampling issues that might come up.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE: Various examples shown...

% runWave(15, 10000, 4, 1, 0, 0, '../data/synthData/ln.roy.10k.wav', '');
% runWave(15, 10000, 3, 1, 1, 0, '../data/synthData/ln.roy.10k.wav', '');
% runWave(15, 10000, 3, 0, 0, 0, '../data/synthData/ln.roy.10k.wav', '../data/synthData/ln.roy.10k.FRM');
% runWave(15, 10000, 3, 1, 0, 0, '../data/synthData/ln.roy.10k.wav', '../data/synthData/ln.roy.10k.FRM');
% runWave(15, 10000, 4, 0, 0, 0, '../data/synthData/ln.roy.10k.wav', '../data/synthData/ln.roy.10k.FRM');
% runWave(15, 10000, 4, 1, 0, 0, '../data/synthData/ln.roy.10k.wav', '../data/synthData/ln.roy.10k.FRM');

% runWave(15, 10000, 3, 0, 1, 1, '../data/synthData/mlm.tea.10k.wav','../data/synthData/mlm.tea.10k.FRM');
% runWave(15, 16000, 4, 1, 1, 0, '../data/VTR_Timit/Timit1.wav', '');
% runWave(15, 16000, 3, 1, 1, 1, '../data/VTR_Timit/Timit1.wav', '../data/VTR_Timit/Timit1.FRM');
% runWave(15, 16000, 4, 1, 1, 1, '../data/VTR_Timit/Timit1.wav', '../data/VTR_Timit/Timit1.FRM');
% runWave(15, 16000, 4, 1, 0, 1, '../data/VTR_Timit/Timit1.wav', '../data/VTR_Timit/Timit1.FRM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../')); %Set paths 
rand('state',sum(100*clock)); randn('state',sum(100*clock)); % Set seeds

%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Checks  %%%%
if(cepOrder < 1)
    error('Must have at least one cepstral coefficient');
end

disp(sprintf('The waveform in %s will be analyzed',dataFileName));
disp(sprintf('The number of formants to be tracked is %d', numFormants));

if(trackBW)
    disp('Tracking bandwidths');
else
    disp('Not tracking bandwidths');
end

if(strcmp(wsFileName,''))
    wvSurf = 0;
    trackBW = 1; % If bandwidths are not provided from wavesurfer, must track them
else
    wvSurf = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER SETTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select which trackers to run, the default is both are set to run
algFlag = [1 1]; % Select 1 to run, 0 not to
EKF = 1; EKS = 2;

% Set sampling frequency
fs = 2*numFormants*1250;

% Voice Activity Detection Parameters
if useVAD
    coastJoint = 1;   % Coast all formants jointly, or not
    quantThresh = .2; % Power quantile threshold for coasting
    plotVad = 0; % Set to '1' for VAD plots
end

% Inter Formant Correlation Parameters
if(wvSurf == 0)
    useCorr = 0; % If no wavesurfer input file is provided, cannot use cross-correlation feature.
end

if useCorr
    diagCorr = 0; % Set to 0 if VAR(1) to 1 when 4 ind AR models
    wsTruth  = 1; % Set to 0 if to estimate matrix from VTR or to 1 for WS
end

% Short-time Analysis Parameters
% These parameters must be MATCHED to the wavesurfer parameters that were
% used to extract formants if such output is provided to the algorithm
wType = 'hamming'; % window type
wLengthMS  = 20;   % Length of window (in milliseconds)
wOverlap = 0.5;    % Factor of overlap of window
lpcOrder = 12;     % Number of LPC Coefficients
peCoeff  = .7;     % Pre-emphasis factor

%%%%%%%%%%%%%%%%%%%%%%%%% Load audio waveform %%%%%%%%%%%%%%%%%%%%%%%%%%%%
wav = wavread(dataFileName);
display(['Reading data from Wavesurfer-generated file ' wsFileName])

%If a wavesurfer file is available read from it
if(wvSurf)
    [trueState, BW_data] = wavesurferFormantRead(wsFileName,numFormants);
end

%%%%%%%%%%%%%%%%%%%%%% Generate observation sequence %%%%%%%%%%%%%%%%%%%%%%

% Resample input data if sampling rate is not equal to that of input
if(fs ~= fs_in)
    wav = resample(wav,fs,fs_in,2048);
    display(['Input Fs = ' int2str(fs_in) ' Hz; resampling to ' num2str(fs) ' Hz'])
end

% Compute window length in samples, now that sampling rate is set
wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);

% Truncate the waveform to match true state
if(wvSurf)
    wav = wav(1:wLength*wOverlap*length(trueState));
end

% Pack the analysis parameters into a datastructure for later use
aParams.wType      = wType;     % Window type
aParams.wLengthMS  = wLengthMS; % Length of window (in milliseconds)
aParams.wLength    = wLength;   % Length of window (in samples)
aParams.win        = win;       % The actual window
aParams.wOverlap   = wOverlap;  % Factor of overlap of window
aParams.lpcOrder   = lpcOrder;  % Number of LPC Coefficients
aParams.peCoeff    = peCoeff;   % Pre-emphasis factor
aParams.fs         = fs;        % Sampling frequency

% Generate cepstral data
y = genLPCC(wav, win, wOverlap, peCoeff, lpcOrder, cepOrder);

% Truncate at most one frame from output of LPCC call
numFrames = size(y,2);
if(wvSurf)
    trueState = trueState(1:numFormants,1:numFrames);
    BW_data = BW_data(1:numFormants,1:numFrames);
end

%%%%%%%%%%%%%%%%%%%%%%% Voice Activity Detection %%%%%%%%%%%%%%%%%%%%%%%%%%
if(useVAD)
    display('Using Voice Activity Detection');
    % If all formants are coasted jointly, then do not do multiband detection
    multiBand = ~coastJoint;
    
    [frameInds] = simpleVAD(wav, aParams, quantThresh, multiBand, [], plotVad);

    if(multiBand)
        if(numFormants > 4)
            error('Multi-band VAD not supported for more than 4 formants');
        else
            % Formant energy bands (taken from Deng et al, 2006)
            bandEdges = [200 900; 600 2800; 1400 3800; 1700 5000];
            bandEdges = bandEdges(1:numFormants,:);
        end

        formantInds = frameInds;
    else
        [frameInds] = simpleVAD(wav, aParams, quantThresh, multiBand, [], plotVad);
        formantInds = repmat(frameInds,numFormants,1)';
    end
else
    display('Not using Voice Activity Detection');
    formantInds = ones(numFrames,numFormants);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Tracking Algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Tracking proceeds via
% x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
% y_{k}   = H_kx_{k} + v_k, v_k ~ N(0, R)

% Need to set the parameters: F, Q and R
% H_k is obtained in the EKF via linearization about the state

% Set process noise covariance matrix Q 
% Currently is estimated from the loaded data from Wavesurfer
Qscale = 1;%.000100; % For experiments scaling this value
if(trackBW)
    if(wvSurf)
        % Get state variance from read in tracks
        Q = diag([var(trueState,0,2); var(BW_data,0,2)])/Qscale;
    else
        % Use hardcoded variance values if wvsurfer not available
        qFormantVar = 500000;
        qBWVar = 5000;
        Q = diag([ones(numFormants,1)*qFormantVar; ones(numFormants,1)*qBWVar]);
    end
else
    Q = diag(var(trueState,0,2))/Qscale;
end

% Set measurement noise covariance matrix R %%%%
% Observation noise variances should decrease as the cepstral order increases
% Using sigExp = 2 for synthetic vowels, but =1 in formant experiments
lambda = 10^(0); sigExp = 1;
R      = diag(1/lambda.*ones(cepOrder,1)./(((1:cepOrder).^sigExp)'));

%%% Process Matrix F %%%%
% Estimate cross-correlation among Formant Trajectories
if(useCorr)
    if(diagCorr)
        for dim = 1:numFormants
            ar(dim) = fitVAR(trueState(dim,:)',1);
        end
        F = diag(ar);
    else
        F = fitVAR(trueState',1)
    end
else
    F = eye(numFormants);
end

if(trackBW)
    F = [F zeros(numFormants); zeros(numFormants) eye(numFormants)];
end

%%% Set Bandwidths %%%%
if(trackBW)
    bwStates = []; % If we are tracking bandwidths do not provide them
else
    bwFlag = 0; % 0 - Use loaded bandwidths, 1 - Average bandwidths
    bwStates = genTrackBW(bwFlag,BW_data);
end

% Initial state of formant trackers
initFormant = 500 + 1000*(0:(numFormants - 1))';
initBW = 80 + 40*(0:(numFormants - 1))';

if(trackBW)
    x0 = [initFormant', initBW']';
else
    x0 = initFormant;
end
display(['Initial state estimates set at: ' num2str(x0')])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now all parameters have been set, go ahead and execute trackers
countTrack = 1; % Counter for storing results

% Initialize root-mean-square error matrices:
rmse    = zeros(numFormants, sum(algFlag));
relRmse = zeros(numFormants, sum(algFlag));

% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;
   [x_estEKF x_errVarEKF] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth);

    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'g-.'};

    countTrack = countTrack + 1;     % Increment counter
end

% Run Extended Kalman Smoother
if algFlag(EKS)
    smooth = 1;  
   [x_estEKS x_errVarEKS] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth);
   
    % Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKS;
    estVar(:,:,:,countTrack) = x_errVarEKS;
    titleCell(1,countTrack+1) = {'EKS'};
    titleCell(2,countTrack+1) = {'b:'};
end

% Basic plotting routines
if(algFlag(EKF))
    % Plot results of Kalman Filter
    plotSpecTracks2(wav, estTracks(1:numFormants,:,1), aParams);
    gcf; title('EKF: Estimated Formant Tracks');

    if(trackBW)
    figure;
    plot(estTracks(numFormants+1:end,:,1)');
    title('Bandwidths Estimated by Extended Kalman Filter (EKF)');
    end
end

if(algFlag(EKS))
    %%%%%%%%%%%%%% Plot results of Kalman Smoother %%%%%%%%%%%%%%%%%%%%%%
    % Overlay plot of formants onto spectrogram
    plotSpecTracks2(wav, estTracks(1:numFormants,:,2), aParams);
    gcf;
    title('EKS: Estimated Formant Tracks');

    % Plot Bandwidths
    if(trackBW)
        figure;
        plot(estTracks(numFormants+1:end,:,2)');
        title('Bandwidths Estimated by Extended Kalman Smoother (EKS)');
    end
end

% If wavesurfer truth is provided, then can compare wavesurfer output to
% estimate tracks
if(wvSurf)
    %Initial Plotting Variables
    titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
    titleCell(2,1)  = {'r'};            % Color for true state plot

    % A basic plotting routine to visualize results, if bandwidths are 
    % tracked their plots will appear with the wavesurfer counterparts
    % after the formant tracks are plotted
    if(trackBW)
        plotStateTracks([trueState; BW_data],estTracks,titleCell);
    else
        plotStateTracks(trueState,estTracks,titleCell);
    end
end