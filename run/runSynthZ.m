function varargout = runSynthZ(F, Fbw, Z, Zbw, N, pNoiseVar, snr_dB, cepOrder, fs, trackBW, plot_flag, algFlag, x0)

% Track poles and zeros (with bandwidths) on synthetic data
% Model-based synthesis, so no speech-like waveform available
% Author: Daniel Rudoy, Daryush Mehta
% Created:  12/13/2007
% Modified: 12/13/2007, 02/15/2010, 03/21/2010
% Modified: 08/23/2010 DDM: track bandwidths, P matrix output
% Modified: 08/25/2010 DDM: pNoiseVar, snr_dB, cepOrder each has two values: one for synthesis, other for tracking
% Modified: 08/31/2010 DDM: aParams deleted as input parameter, trueState
%                           added as output, observation noise variance for tracker is
%                           calculated from snr_dB as if we had total
%                           number of cepstral coefficients, cepOrder(2)
%                           then truncates that vector
% 
% INPUT:
%    F:         center frequencies of the resonances (col vector), in Hz
%    Fbw:       corresponding bandwidths of the resonances (col vector), in Hz
%    Z:         center frequencies of the anti-resonances (col vector), in Hz
%    Zbw:       corresponding bandwidths of the anti-resonances (col vector), in Hz
%    N:         number of observations to generate
%    pNoiseVar: process noise variance; first element for synthesis, second
%                     element for tracking; [100 120]
%    snr_dB:    observation noise, in dB; first element for synthesis, second
%                     element for tracking; [15 25]
%    cepOrder:  Number of cepstal coefficients to compute; first element for synthesis, second
%                     element for tracking; [15 25]
%    fs:        sampling rate of waveform, in Hz [no waveform generated but fs used in
%                     mapping from formant parameters to ARMA cepstrum]
%    trackBW:   track bandwidths if 1
%    plot_flag: plot figures if 1
%    algFlag:   select 1 to run, 0 not to for [EKF EKS]
%    x0:        initial state of formant trackers [F; Fbw; Z; Zbw], in Hz
% 
% OUTPUT:
%    Depends on algFlag. For each algorithm, three outputs plus the trueState--
%       rmse:       RMSE for each track
%       x_est:      estimated tracks
%       x_errVar:   covariance matrix of estimated tracks
%       So that if two algorithms run, the following are output:
%       [rmseEKF, x_estEKF, x_errVarEKF, rmseEKS, x_estEKS, x_errVarEKS, trueState]
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% see runSynthZ_wrapper
%
% The usual examples from PRAAT and Wavesurfer are not included, because
% it is not clear how to generate paths of zeros along with the paths of
% formants that are available.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath(genpath('../')); % Paths

%% %%%%%%%%%%%%%%%%%%% Synthesis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initState    = [F; Z];
initBW       = [Fbw; Zbw]';

%% Generate data (this simply evolves center frequencies and bandwidths via random walk)
[trueState BW_data] = genSynthFormantTracks(pNoiseVar(1), N, initState, initBW);

%% Generate noisy LPCC coefficients (observation sequence)
[y_synth, oNoiseVar, C_clean] = genNoisyObserZ(snr_dB(1), trueState(1:length(F),:), BW_data(1:length(F),:), trueState(length(F)+1:end,:), BW_data(length(F)+1:end,:), cepOrder(1), fs);

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
y = y_synth(1:cepOrder(2), :);
numObs = size(y, 2);

if trackBW
    trueState = [trueState(1:length(F),:); BW_data(1:length(F),:); trueState(length(F)+1:end,:); BW_data(length(F)+1:end,:)];
    bwStates = []; % If we are tracking bandwidths do not provide them
    Qdiag = pNoiseVar(2)*[ones(1,length(F)), 1/25*ones(1,length(F)), ones(1,length(Z)), 1/25*ones(1,length(Z))];
else
    bwStates = BW_data;
    Qdiag = pNoiseVar(2)*ones(1,length(F)+length(Z));
end
stateLen = size(trueState,1);

% Process Matrix F
Fmatrix = eye(stateLen);

% process noise from read in tracks
Q = diag(Qdiag); % use input variance for model
% Q = diag(var(trueState,0,2));
% Q = diag(pNoiseVar(2)*ones(1, stateLen));
% Q = diag([1000 1000 1000 1000]);

randn('state',sum(100*clock)); % different start point each time the noise is created
% randn('state',5); % so added noise is the same each time
[Ctemp, noiseVar] = addONoise(C_clean(1:cepOrder(1), :), snr_dB(2)); clear Ctemp
R = noiseVar*eye(cepOrder(2)); % Measurement noise covariance matrix R

% A voice activity detector is not used here in the synthetic case
formantInds = ones(stateLen,numObs);

countTrack = 1; % Counter for storing results
countOut = 1; % Counter for output variables

EKF = 1; EKS = 2;

% Initialize root-mean-square error matrices
rmse    = zeros(stateLen, sum(algFlag));
relRmse = zeros(stateLen, sum(algFlag));

% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;
    [x_estEKF x_errVarEKF] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, length(F), smooth);

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
    [x_estEKS x_errVarEKS] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, length(F), smooth);

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

varargout(countOut) = {trueState}; countOut = countOut + 1;

if plot_flag
    %% Initial Plotting Variables
    titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
    titleCell(2,1)  = {'r'};            % Color for true state plot

    % A basic plotting routine to visualize results
    for j = 1:sum(algFlag)
        plotStateTracksFZ(trueState,estTracks(:,:,j),titleCell(:,[1 j+1]), length(F),trackBW);
    end
end