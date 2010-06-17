function varargout = runWaveZ(cepOrder, numFormants, numAntiF, trackBW, dataFileName, algFlag, x0)

% This code is incomplete. It is meant to test the pole/zero trackers on
% real speech. However, two things are necessary
% a) Bandwidths supplied for formants (this is clear how to do)
% b) Bandwidths supplied for zeros    (this is less clear how to do)
% One possibility would be to use the same bw values for zeros as for
% formants. That seems 'reasonable'. Outside of that, we will require to
% write a tracker for zeros and bandwidths simultaneously. This has not
% yet been done.
% DGR: 02/15/2010

% UPDATE
% Tracker has been coded for zeros and bandwidths simultaneously.
% DDM: 06/02/2010

% Test tracking algorithms using raw audio waveform. Requires a 
% post-processed file from WaveSurfer in order to input bandwidths,
% estimate inter-formant correlation (if desired) and validate the results.
% Nothing here depends on the VTR database, if you wish to use the VTR
% database for testing use testVTR.m
%
% Author: Daniel Rudoy, Daryush Mehta
% Created:  12/13/2007
% Modified: 01/07/2007, 06/02/10

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
%    x0:        initial states (center frequencies) of formant trackers [(2*numFormants + 2*numAntiformants) x 1], in Hz
% 
% OUTPUT:
%    Depends on algFlag. For each algorithm, two outputs generated--
%       rmse_mean1: average RMSE across all tracks
%       x_est1:  estimated tracks
%       So that if two algorithms run, the following are output:
%       [rmse_mean1, x_est1, rmse_mean2, x_est2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% see runWaveZ_wrapper.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Error checking

% Set paths
addpath(genpath('../'));

% Set random seeds to allow for exact repetition or otherwise
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

%%% Load audio waveform %%%
display(['Reading data from WAV file ' dataFileName])
[wav, fs_in] = wavread(dataFileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Generate observation sequence %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analysis parameters
wType = 'hanning'; % Window type
wLengthMS  = 100;   % Length of window (in milliseconds)
wOverlap = 0.5;    % Factor of overlap of window
lpcOrder = 10;     % Number of LPC Coefficients
zOrder = 4;        % Number of MA coefficients
peCoeff  = .9;     % Pre-emphasis factor

% Set sampling rate
fs = 3000;

% Resample input data if sampling rate is not equal to that of input
if(fs ~= fs_in)
    display(['Input Fs = ' int2str(fs_in) ' Hz; resampling to ' num2str(fs) ' Hz'])
    wav = resample(wav,fs,fs_in,2048);
end

% Compute window length in samples, now that sampling rate is set
wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);

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
y = genLPCCz(wav, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Tracking Algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for tracking algorithms %%%
% x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
% y_{k}   = Hx_{k} + v_k, v_k ~ N(0, R)
% We need to set the parameters: F, Q and R
% H is obtained in the EKF via linearization about the state
nP = numFormants;
numObs = size(y,2);

if trackBW
    stateLen = 2*numFormants + 2*numAntiF;
    bwStates = []; % If we are tracking bandwidths do not provide them
else
    stateLen = numFormants + numAntiF;
    
    bwStates = repmat([40.56219      7.451091   50]', 1, numObs);
end

% A voice activity detector is not used here yet
display('Not using Voice Activity Detection');
formantInds = ones(stateLen, numObs);

% Process Matrix F, Correlation is not being tested here
Fmatrix = eye(stateLen);

% Process noise covariance matrix Q
Q = 2500*eye(stateLen); %%% NEED BETTER ESTIMATION

% Measurement noise covariance matrix R
% Set/choose estimated observation noise variance, which should decrease
% as the cepstral order increases
% Using sigExp = 2 for synthetic vowels, but =1 in formant experiments
lambda = 10^(0); sigExp = 1;
R = diag(1/lambda.*ones(cepOrder,1)./(((1:cepOrder).^sigExp)'));

countTrack = 1; % Counter for storing results
countOut = 1; % Counter for output variables

% Select which trackers to run
EKF = 1; EKS = 2;

% Initialize root-mean-square error matrices:
rmse    = zeros(numFormants, sum(algFlag));
relRmse = zeros(numFormants, sum(algFlag));

%% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 1;
    [x_estEKF x_errVarEKF] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);
    
    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    
    varargout(countOut) = {x_estEKF}; countOut = countOut + 1;
    varargout(countOut) = {x_errVarEKF}; countOut = countOut + 1;
end

% Run Extended Kalman Smoother
if algFlag(EKS)
    smooth = 1;
    [x_estEKS x_errVarEKS] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);    
    
    % Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKS;
    estVar(:,:,:,countTrack) = x_errVarEKS;
    
    varargout(countOut) = {x_estEKS}; countOut = countOut + 1;
    varargout(countOut) = {x_errVarEKS}; countOut = countOut + 1;
end

varargout(countOut) = {aParams};

% Super-impose over a spectrogram
% plotSpecTracks2(wav, estTracks(1:numFormants+numAntiF,:), aParams, numAntiF);