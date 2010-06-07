function [] = runWaveZ(cepOrder, fs_in, numFormants, numZeros, trackBW, dataFileName, wsFileName)

% This code is incomplete. It is meant to test the pole/zero trackers on
% real speech. However, two things are necessary
% a) Bandwidths supplied for formants (this is clear how to do)
% b) Bandwidths supplied for zeros    (this is less clear how to do)
% One possibility would be to use the same bw values for zeros as for
% formants. That seems 'reasonable'. Outside of that, we will require to
% write a tracker for zeros and bandwidths simultaneously. This has not
% yet been done.
% DGR: 02/15/2010

% Test tracking algorithms using raw audio waveform. Requires a 
% post-processed file from WaveSurfer in order to input bandwidths,
% estimate inter-formant correlation (if desired) and validate the results.
% Nothing here depends on the VTR database, if you wish to use the VTR
% database for testing use testVTR.m
%
% Author: Daniel Rudoy
% Created:  12/13/2007
% Modified: 01/07/2007

% INPUT
%
% cep_order :   How many cepstal coefficients to include in the observations
% fs_in     :   Sampling rate at which the observations are made (fake here)
% numFormants: Number of formants that we should track
% dataFileName: Audio (.wav) source file
% wsFileName  : Currently compared against Wavesurfer. This must contain
%               at least as many formants as we wish to track and
%               preferrably the exact number to avoid mismatch in
%               resampling issues that might come up.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% runWaveZ(15, 10000, 4, 0, 0,  '../data/synthData/mlm.tea.10k.wav', '../data/synthData/mlm.tea.10k.FRM');
% runWaveZ(15, 10000, 3, 0, 0, '../data/synthData/ln.roy.10k.wav', '../data/synthData/ln.roy.10k.FRM');
% runWaveZ(15, 10000, 4, 0, 0, '../data/synthData/ln.roy.10k.wav', '../data/synthData/ln.roy.10k.FRM');
% runWaveZ(15, 16000, 3, 0, 0, '../data/VTR_Timit/Timit1.wav', '../data/VTR_Timit/Timit1.FRM');
% runWaveZ(15, 10000, 4, 0, 1, '../data/synthData/mlm.tea.10k.wav', '../data/synthData/mlm.tea.10k.FRM');
% runWaveZ(15, 10000, 3, 0, 1, '../data/synthData/ln.roy.10k.wav', '../data/synthData/ln.roy.10k.FRM');
% runWaveZ(15, 10000, 4, 0, 1, '../data/synthData/ln.roy.10k.wav', '../data/synthData/ln.roy.10k.FRM');
% runWaveZ(15, 16000, 3, 0, 1, '../data/VTR_Timit/Timit1.wav', '../data/VTR_Timit/Timit1.FRM');
% runWaveZ(15, 16000, 3, 0, 1, '../data/synthData/female_e.wav', '../data/synthData/female_e.frm')
% runWaveZ(15, 16000, 3, 0, 1, '../data/synthData/female_e_noise.wav', '../data/synthData/female_e_noise.frm')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Error checking

% Set paths
addpath(genpath('../'));

% Set random seeds to allow for exact repetition or otherwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('rand_state','var')
    rand_state = sum(100*clock);
    display('Seeding rand generator according to clock')
else
    display('Resetting rand generator according to fixed input')
end
if ~exist('randn_state','var')
    pause(rand); % Try to decorrelate rand seed from randn seed
    randn_state = sum(100*clock);
    display('Seeding randn generator according to clock')
else
    display('Resetting randn generator according to fixed input')
end
rand('state',rand_state);
randn('state',randn_state);


% Correlation control
useCorr = 0; % Flag for whether or not to use correlation
if useCorr
    diagCorr = 0; % Set to 0 if VAR(1) to 1 when 4 ind AR models
    wsTruth  = 1; % Set to 0 if to estimate matrix from VTR or to 1 for WS
end

% Flag for controlling whether or not to use voice activity detection
% to coast the formant tracks or for some other purpose
useVAD = 0;
if useVAD
    coastJoint = 1;   % Coast all formants jointly, or not
    quantThresh = .2; % power quantile threshold for coasting
end


%%% Load audio waveform %%%
wav = wavread(dataFileName);
display(['Reading data from Wavesurfer-generated file ' wsFileName])
[trueState, BW_data] = wavesurferFormantRead(wsFileName,numFormants);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Generate observation sequence %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analysis parameters
wType = 'hamming'; % Window type
wLengthMS  = 20;   % Length of window (in milliseconds)
wOverlap = 0.5;    % Factor of overlap of window
lpcOrder = 12;     % Number of LPC Coefficients
peCoeff  = .7;     % Pre-emphasis factor

% Set sampling rate
% This yields 3750/5000 Hz if using 3/4 formants respectively
fs = 2*numFormants*1250;
%fs = 2*3750;

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
wav = wav(1:wLength*wOverlap*length(trueState));

% Pack the analysis parameters for later use
aParams.wType      = wType;     % Window type
aParams.wLengthMS  = wLengthMS; % Length of window (in milliseconds)
aParams.wLength    = wLength;   % Length of window (in samples)
aParams.win        = win;       % The actual window
aParams.wOverlap   = wOverlap;  % Factor of overlap of window
aParams.lpcOrder   = lpcOrder;  % Number of LPC Coefficients
aParams.peCoeff    = peCoeff;   % Pre-emphasis factor
aParams.fs         = fs;        % Sampling frequency

% Generate cepstral data
% No longer require the function genLPCC 2 since VAD has been pulled out
y = genLPCC(wav, win, wOverlap, peCoeff, lpcOrder, cepOrder);

% for i = 1:15
%     figure(i);
%     plot(y(i,:));
% end
% Truncate at most one frame from output of LPCC call
numFrames = size(y,2);
numObs = numFrames;
trueState = trueState(1:numFormants,1:numFrames);
BW_data = BW_data(1:numFormants,1:numFrames);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Voice Activity Detection %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(useVAD)
    display('Using Voice Activity Detection');
    % If all formants are coasted jointly, then do not do multiband detection
    multiBand = ~coastJoint;
    % Set some parameters prior to calling the VAD
    plotVad = 1;

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

%%% Set parameters for tracking algorithms %%%

% x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
% y_{k}   = Hx_{k} + v_k, v_k ~ N(0, R)

% We need to set the parameters: F, Q and R
% H is obtained in the EKF via linearization about the state

%%% Process noise covariance matrix Q %%%%

% Currently is estimated from the loaded data from Wavesurfer
Qscale = 1;%.000100; % For experiments scaling this value
if(trackBW)
    Q = diag([var(trueState,0,2); var(BW_data,0,2)])/Qscale;     % Get state variance from read in tracks    
else
    Q = diag(var(trueState,0,2))/Qscale;
end

%%% Measurement noise covariance matrix R %%%%

% Set/choose estimated observation noise variance, which should decrease
% as the cepstral order increases
% Using sigExp = 2 for synthetic vowels, but =1 in formant experiments
lambda = 10^(0); sigExp = 1;
R      = diag(1/lambda.*ones(cepOrder,1)./(((1:cepOrder).^sigExp)'));

%%% Process Matrix F %%%%

% Correlation is not being tested here

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

% The formant trackers here do not themselves track bandwidths, so an
% appropriate array of bandwidths based on say Wavesurfer input is required
bwFlag = 0; % 0 - Use loaded bandwidths, 1 - Average bandwidths
bwStates = genTrackBW(bwFlag,BW_data);

% Initial state of formant trackers
initFormant = 500 + 1000*(0:(numFormants - 1))';
%initFormant = [500 500 1500]';
initBW = 80 + 40*(0:(numFormants - 1))';

if(trackBW)
    x0 = [initFormant', initBW']';
else
    x0 = initFormant; 
end
display(['Initial state estimates set at: ' num2str(x0')])

% Select which trackers to run
algFlag = [1 1]; % Select 1 to run, 0 not to
EKF = 1; EKS = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now all parameters have been set, go ahead and execute trackers

countTrack = 1; % Counter for storing results

% Initialize root-mean-square error matrices:
rmse    = zeros(numFormants, sum(algFlag));
relRmse = zeros(numFormants, sum(algFlag));

% Run Extended Kalman Filter
if algFlag(EKF)

    t0 = clock;     % Run and time EK filter

    if(trackBW)
        [x_estEKF x_errVarEKF] = formantTrackEKFBW(y, F, Q, R, x0, formantInds, fs);
    else
        [x_estEKF x_errVarEKF] = formantTrackEKF(y, F, Q, R, x0, formantInds, fs, bwStates);
    end
    
    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'g-.'};

    % Compute and Display MSE and RMSE
    for j = 1:numFormants
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
    end

    % Display output summary and timing information
    display(['Kalman Filter Run Time: ' num2str(etime(clock,t0)) ' s. Average RMSE: ' num2str(mean(rmse(:,countTrack)))]);
    countTrack = countTrack + 1;     % Increment counter
end

% Run Extended Kalman Smoother
if algFlag(EKS)

    t0 = clock;
 
    if(trackBW)
        [x_estEKS x_errVarEKS] = formantTrackEKSBW(y, F, Q, R, x0, formantInds, fs);
    else
        [x_estEKS x_errVarEKS] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates);
    end
    
    
    % Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKS;
    estVar(:,:,:,countTrack) = x_errVarEKS;
    titleCell(1,countTrack+1) = {'EKS'};
    titleCell(2,countTrack+1) = {'b:'};

    %Compute and Display MSE and RMSE
    for j = 1:numFormants
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
    end

    % Display output summary and timing information
    display(['Kalman Smoother Run Time: ' num2str(etime(clock,t0)) ' s. Average RMSE: ' num2str(mean(rmse(:,countTrack)))]);
    countTrack = countTrack + 1;    % Increment counter
end


%Initial Plotting Variables
titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};            % Color for true state plot

% A basic plotting routine to visualize results
if(trackBW)
    plotStateTracks([trueState; BW_data],estTracks,titleCell);
else
    plotStateTracks(trueState,estTracks,titleCell);
end

% Super-impose over a spectrogram
plotSpecTracks2(wav, estTracks(1:numFormants,:), aParams);
plotSpecTracks2(wav, trueState, aParams);