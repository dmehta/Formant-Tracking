% specifically called by runWaveZ_wrapper_nan_synth9v.m

function varargout = runWaveZ(cepOrder, numFormants, numAntiF, trackBW, dataFileName, algFlag, x0, bwStates, aParams, formantInds)

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
% DDM: 06/02/2010 Tracker has been coded for zeros and bandwidths simultaneously.
% DDM: 09/28/2010 formantInds is input parameter
% 
% Test tracking algorithms using raw audio waveform. Requires a 
% post-processed file from WaveSurfer in order to input bandwidths,
% estimate inter-formant correlation (if desired) and validate the results.
% Nothing here depends on the VTR database, if you wish to use the VTR
% database for testing use testVTR.m
%
% Author: Daniel Rudoy, Daryush Mehta
% Created:  12/13/2007
% Modified: 01/07/2007, 06/02/10
%           06/16/2010 (analysis parameters on input)
% 
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
%    x0:            Initial states (center frequencies) of formant trackers [(2*numFormants + 2*numAntiformants) x 1], in Hz
%    bwStates:      Initial states (bandwidths) of formant trackers, if
%                       applicable [(2*numFormants + 2*numAntiformants) x 1], in Hz
%    aParams:       Structure with following fields:
%                       wType       % window type
%                       wLengthMS   % Length of window (in milliseconds)
%                       wOverlap    % Factor of overlap of window
%                       lpcOrder    % Number of AR coefficients
%                       zOrder      % Number of MA coefficients
%                       peCoeff     % Pre-emphasis factor
%    formantInds:   masking matrix of formant indices; 0 if masked, 1 if unmasked
% 
% OUTPUT:
%    Depends on algFlag. For each algorithm, two outputs then waveform--
%       rmse_mean:  average RMSE across all tracks
%       x_est:      estimated tracks
%       x_errVar:   covariance matrix of estimated tracks
%       So that if two algorithms run, the following are output:
%       [rmse_meanEKF, x_estEKF, x_errVarEKF, rmse_meanEKS, x_estEKS, x_errVarEKS, wav]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% see runWaveZ_wrapper.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Error checking

% Set paths
% addpath(genpath('../'));

% Load audio waveform %
display(['Reading data from WAV file ' dataFileName])
[wav, fs_in] = wavread(dataFileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Generate observation sequence %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analysis parameters
wType       = aParams.wType;     % Window type
wLengthMS   = aParams.wLengthMS; % Length of window (in milliseconds)
wOverlap    = aParams.wOverlap;  % Factor of overlap of window
lpcOrder    = aParams.lpcOrder;  % Number of AR Coefficients
zOrder      = aParams.zOrder;    % Number of MA Coefficients
peCoeff     = aParams.peCoeff;   % Pre-emphasis factor
fs          = aParams.fs;        % Sampling frequency

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
aParams.wLength    = wLength;   % Length of window (in samples)
aParams.win        = win;       % The actual window

% Generate cepstral data
y = genLPCCz(wav, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);
% y = genCeps(wav, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);

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
    numAntiF = length(x0)-numFormants;
    x0 = [x0(1:numFormants); bwStates(1:numFormants); ...
        x0(numFormants+1:numFormants+numAntiF); bwStates(numFormants+1:numFormants+numAntiF)];
    bwStates = []; % If we are tracking bandwidths do not provide them

    % Process noise covariance matrix Q
    Q = zeros(size(x0));
    %Q = 100^2*eye(stateLen); %%% NEED BETTER ESTIMATION
    Q(1:numFormants) = 100^2;
    Q(numFormants+1:numFormants*2) = 50^2;
    Q(2*numFormants+1:2*numFormants+numAntiF) = 100^2;
    Q(2*numFormants+numAntiF+1:end) = 50^2;
    Q = diag(Q);
else
    stateLen = numFormants + numAntiF;
    bwStates = repmat(bwStates, 1, numObs);
    
    % Process noise covariance matrix Q
    Q = diag(ones(stateLen,1)*10e4);
end

% formantInds is input
if isempty(formantInds)
    display('Not using Voice Activity Detection');
    formantInds = ones(stateLen, numObs);
else
%     coastJoint = 1;    % Coast all formants jointly, or not
%     quantThresh = .15; % power quantile threshold for coasting
%     plotVad = 0;
%     multiBand = ~coastJoint;
%     [frameInds] = simpleVAD(wav, aParams, quantThresh, multiBand, [], plotVad);
%     if trackBW
%         formantInds = repmat(frameInds,(numAntiF+2)*numFormants,1)';
%     else
%         formantInds = repmat(frameInds,numFormants+numAntiF,1)';
%     end
%     formantInds = formantInds';
end

% Process Matrix F, Correlation is not being tested here
Fmatrix = eye(stateLen);

% Measurement noise covariance matrix R
% Set/choose estimated observation noise variance, which should decrease
% as the cepstral order increases
% Using sigExp = 2 for synthetic vowels, but =1 in formant experiments
lambda = 10^(0); sigExp = 1;
R = diag(1/lambda.*ones(cepOrder,1)./(((1:cepOrder).^sigExp)'));
% R = diag(1e-0*ones(cepOrder,1));

countTrack = 1; % Counter for storing results
countOut = 1; % Counter for output variables

% Select which trackers to run
EKF = 1; EKS = 2;

% Initialize root-mean-square error matrices:
rmse    = zeros(numFormants, sum(algFlag));
relRmse = zeros(numFormants, sum(algFlag));

%% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;
    [x_estEKF x_errVarEKF] = formantTrackEKSZ_nan_synth9v(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);
    
    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    
    varargout(countOut) = {x_estEKF}; countOut = countOut + 1;
    varargout(countOut) = {x_errVarEKF}; countOut = countOut + 1;
end

% Run Extended Kalman Smoother
if algFlag(EKS)
    smooth = 1;
    [x_estEKS x_errVarEKS] = formantTrackEKSZ_nan_synth9v(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);    
    
    % Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKS;
    estVar(:,:,:,countTrack) = x_errVarEKS;
    
    varargout(countOut) = {x_estEKS}; countOut = countOut + 1;
    varargout(countOut) = {x_errVarEKS}; countOut = countOut + 1;
end

varargout(countOut) = {wav}; countOut = countOut + 1;

% Super-impose over a spectrogram
% plotSpecTracks2(wav, estTracks(1:numFormants+numAntiF,:), aParams, numAntiF);