function varargout = runSynth_OLA(F, Fbw, Z, Zbw, N, snr_dB, cepOrder, fs, trackBW, plot_flag, algFlag, x0, aParams, sParams, f0info)

% Track formants and anti-formants (no bandwidths) on synthetic data
% that constructs a waveform using overlap-add of windows generated from an
% ARMA model.
% 
% Author: Daryush Mehta
% Created:  05/09/2010
% Modified: 05/17/2010 (bandwidth contours)
%           06/16/2010 (analysis parameters on input, x output)
%           08/31/2010 DDM: snr_dB, cepOrder each has two values: one for synthesis, other for tracking
%           12/13/2010 DDM: f0 as input and voicing indicator and speech
%           indicator
% 
% INPUT:
%    F:         center frequencies of the resonances (numFormants x numFrames), in Hz
%    Fbw:       corresponding bandwidths of the resonances (numFormants x numFrames), in Hz
%    Z:         center frequencies of the anti-resonances (numAntiformants x numFrames), in Hz
%    Zbw:       corresponding bandwidths of the anti-resonances (numAntiformants x numFrames), in Hz
%    N:         length of signal, in samples
%    snr_dB:    observation noise, in dB; first element for synthesis, second
%                     element for tracking; [15 25]
%    cepOrder:  Number of cepstal coefficients to compute; first element for synthesis, second
%                     element for tracking; [15 25]
%    fs:        sampling rate of waveform, in Hz
%    trackBW:   track bandwidths if 1
%    plot_flag: plot figures if 1
%    algFlag:   select 1 to run, 0 not to for [EKF EKS]
%    x0:        initial states (center frequencies) of formant trackers [(2*numFormants + 2*numAntiformants) x 1], in Hz
%    aParams:   structure with following fields:
%                     wType       % window type
%                     wLengthMS   % Length of window (in milliseconds)
%                     wOverlap    % Factor of overlap of window
%                     lpcOrder    % Number of AR coefficients
%                     zOrder      % Number of MA coefficients
%                     peCoeff     % Pre-emphasis factor
%    sParams:   structure with following fields:
%                     wType       % window type
%                     wLengthMS   % Length of window (in milliseconds)
%                     wOverlap    % Factor of overlap of window
%    f0info:    speech/silence (1/0), voicing (1/0), and f0 (Hz) for each frame (3 x numFrames)
% 
% OUTPUT:
%    Depends on algFlag. For each algorithm, three outputs then waveform, then true state--
%       rmse:       RMSE for each track
%       x_est:      estimated tracks
%       x_errVar:   covariance matrix of estimated tracks
%       trueState:  true state matrix
%       So that if two algorithms run, the following are output:
%       [rmseEKF, x_estEKF, x_errVarEKF, rmseEKS, x_estEKS, x_errVarEKS, x, trueState]
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% Synthetic Examples:
%   see runSynth_OLA_wrapper.m

% addpath(genpath('../')); % Paths
rand('state',sum(100*clock)); randn('state',sum(100*clock)); % Seeds
% rand('state', 2); randn('state', 44);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(N, 1);
windows = zeros(N, 1);
clear allCoeffsP allCoeffsZ

wType       = sParams.wType;     % Window type
wLengthMS   = sParams.wLengthMS; % Length of window (in milliseconds)
wOverlap    = sParams.wOverlap;  % Factor of overlap of window

wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);

% Compute left and right boundaries
wLeft  = 1:wLength*(1-wOverlap):N-wLength+1;
wRight = wLength:wLength*(1-wOverlap):N;

% Compute number of frames
if ~isempty(F)
    numFrames = size(F, 2);
else
    numFrames = size(Z, 2);
end

% If input is scalar, expand to number of frames
if size(F, 2) == 1
    F = repmat(F, 1, numFrames);
end
if size(Z, 2) == 1
    Z = repmat(Z, 1, numFrames);
end

if size(Fbw, 2) == 1
    Fbw = repmat(Fbw, 1, numFrames);
end
if size(Zbw, 2) == 1
    Zbw = repmat(Zbw, 1, numFrames);
end

% % output true state
% if trackBW
%     trueState = [F; Fbw; Z; Zbw];
% else
%     trueState = [F; Z]; % 
% end

%%
for i=1:numFrames
    if ~isempty(F) && ~isempty(Z)
        [num, denom] = fb2tf(F(:, i), Fbw(:, i), Z(:, i), Zbw(:, i), fs);
    end
    
    if ~isempty(F) && isempty(Z)
        [num, denom] = fb2tf(F(:, i), Fbw(:, i), [], [], fs);
    end
    
    if isempty(F) && ~isempty(Z)
        [num, denom] = fb2tf([], [], Z(:, i), Zbw(:, i), fs);
    end

    sInd = f0info(1, i);
    vInd = f0info(2, i);
    f0   = f0info(3, i);
    if sInd == 0 % silence frame
        tmp = 0.1*randn(wLength,1); %% ZERO OR SMALL NOISE??
    else
        if vInd > 0 % then voiced
            if f0 == 0 % catch error
                f0 = 100; % default
            end
            noise_type = '1';
            vowel = 'a';
            gender = 'm';
            dur = wLengthMS/1000;
            dcFlow = 0;
            OQ = 60;
            jitter = 0;
            shimmer = 0;
            HNR = 20;
            timeOffset = 0;
            formants_bw = [80 80 80];
            formants_freq = [730 1090 2440];
            source_type = '3'; % 3 = derivative
            tilt = 0;
            SQ = 270;
            [output, fs, periodicSynth, noiseSynth, inputSource, noiseSource, actualPitchContour, ...
                periodicSynthNoRadiation, noiseSynthNoRadiation, ARcoefs] = ...
                vowelSynth_tilt_arbitrarysource_sourceHNR(noise_type, vowel, f0, gender, fs, dur, dcFlow, OQ, ...
                jitter, shimmer, HNR, timeOffset, formants_freq, formants_bw, source_type, tilt, SQ);

            % choose one:
            inputSource = inputSource';
            %inputSource = diff(inputSource)';
            % inputSource = noiseSource';

            tmp = filter(num, denom, inputSource); % periodic input
        else % unvoiced
            tmp = filter(num, denom, randn(wLength,1)); % stochastic input
        end
    end

    x(wLeft(i):wRight(i)) = x(wLeft(i):wRight(i)) + win.*tmp;
    windows(wLeft(i):wRight(i)) = windows(wLeft(i):wRight(i)) + win;

    allCoeffsP(:, i) = denom;
    allCoeffsZ(:, i) = num;
end

x = (x-mean(x));

% Generate noisy LPCC coefficients (observation sequence), y_synth not used
% but available to compare with C_clean to see impact of adding noise in
% waveform domain vs noise in cepstral domain
[y_synth, oNoiseVar, C_clean] = genNoisyObserZ(snr_dB(1), F, Fbw, Z, Zbw, cepOrder(1), fs);

% NB: could also, equivently, obtain cepstral coefficients using ARMA coefficients
% C_clean2 = lpc2cz(-allCoeffsP(2:end,:),-allCoeffsZ(2:end,:),cepOrder(1));
% sum(sum(C_clean-C_clean2)); % should be zero

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate LPCC coefficients in each window
wType       = aParams.wType;     % Window type
wLengthMS   = aParams.wLengthMS; % Length of window (in milliseconds)
wOverlap    = aParams.wOverlap;  % Factor of overlap of window
lpcOrder    = aParams.lpcOrder;  % Number of AR Coefficients
zOrder      = aParams.zOrder;    % Number of MA Coefficients
peCoeff     = aParams.peCoeff;   % Pre-emphasis factor

wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);

y = genLPCCz(x, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder(2)); % ARMA cepstrum
% y = genCeps(x, win, wOverlap, peCoeff, cepOrder(2)); % cepstrum

%% Do a plot of the LPCCC observations
if plot_flag
    figure;
    imagesc(log(abs(y))); colorbar;
    title('Cepstral Coefficients');
    xlabel('Frame Number');
end

%% %%%%%%%%%%%%%%%%%%% Tracking Algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for tracking algorithms %%%
% x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
% y_{k}   = Hx_{k} + v_k, v_k ~ N(0, R)
% We need to set the parameters: F, Q and R
% H is obtained in the EKF via linearization about the state
nP = size(F, 1);
numObs = size(y,2);

if trackBW
    trueState = [F; Fbw; Z; Zbw];
    bwStates = []; % If we are tracking bandwidths do not provide them
else
    trueState = [F; Z];
    bwStates = [Fbw; Zbw];
end
stateLen = size(trueState,1);

% Process Matrix F
Fmatrix = eye(stateLen);

% process noise estimated from variance of known data
Q = diag(var(trueState(:,2:end)-trueState(:,1:end-1),0,2)+eps);
% Q = Q.*diag([1 1 100 100 1 100]);

[Ctemp, noiseVar] = addONoise(C_clean(1:cepOrder(2), :), snr_dB(2)); clear Ctemp
R = noiseVar*eye(cepOrder(2)); % Measurement noise covariance matrix R

% A voice activity detector is not used here in the synthetic case
formantInds = ones(stateLen, numObs);

countTrack = 1; % Counter for storing results
countOut = 1; % Counter for output variables

EKF = 1; EKS = 2;

% Initialize root-mean-square error matrices
rmse    = zeros(stateLen, sum(algFlag));
relRmse = zeros(stateLen, sum(algFlag));

% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;
    [x_estEKF x_errVarEKF] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);

    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'b:'};

    % Compute and Display MSE and RMSE
    for j = 1:stateLen
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
    end

    % Display output summary and timing information
    rmse_mean = mean(rmse(:,countTrack));
    varargout(countOut) = {rmse}; countOut = countOut + 1;
    varargout(countOut) = {x_estEKF}; countOut = countOut + 1;
    varargout(countOut) = {x_errVarEKF}; countOut = countOut + 1;
    
    countTrack = countTrack + 1;     % Increment counter
end

% Run Extended Kalman Smoother
if algFlag(EKS)
    smooth = 1;
    [x_estEKS x_errVarEKS] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);

    % Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKS;
    estVar(:,:,:,countTrack) = x_errVarEKS;
    titleCell(1,countTrack+1) = {'EKS'};
    titleCell(2,countTrack+1) = {'b:'};

    % Compute and Display MSE and RMSE
    for j = 1:stateLen
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
    end

    % Display output summary and timing information
    rmse_mean = mean(rmse(:,countTrack));
    varargout(countOut) = {rmse}; countOut = countOut + 1;
    varargout(countOut) = {x_estEKS}; countOut = countOut + 1;
    varargout(countOut) = {x_errVarEKS}; countOut = countOut + 1;
    
    countTrack = countTrack + 1;     % Increment counter
end

varargout(countOut) = {x}; countOut = countOut + 1;
varargout(countOut) = {trueState}; countOut = countOut + 1;

if plot_flag
    %% Initial Plotting Variables
    titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
    titleCell(2,1)  = {'r'};            % Color for true state plot

    % A basic plotting routine to visualize results
    for j = 1:sum(algFlag)
        plotStateTracksFZ(trueState,estTracks(:,:,j),titleCell(:,[1 j+1]), nP, trackBW);
    end
end