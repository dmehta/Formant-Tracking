function rmse = runSynthZ(nP, nZ, N, pNoiseVar, snr_dB, cepOrder, fs)

% Track poles and zeros (no bandwidths) and on synthetic data
% Author: Daniel Rudoy
% Created:  12/13/2007
% Modified: 12/13/2007, 02/15/2010, 03/21/2010
% 
% INPUT
%
% nP  : Number of pole pairs to track
% nZ  : Number of zero pairs to track
% N    : Number of observations to generate
% pNoiseVar : Process noise variance
% snr_dB    : How much observation noise to add
% cepOrder : How many cepstal coefficients to include in the observations
% fs        : Sampling rate at which the observations are made (fake here)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% Synthetic Examples:
%   runSynthZ(4,4,100,30^2, 15, 15, 16000);
%   runSynthZ(4,2,250,20^2, 10, 15, 10000);
%
% The usual examples from PRAAT and Wavesurfer are not included, because
% it is not clear how to generate paths of zeros along with the paths of
% formants that are available.

addpath(genpath('../')); % Paths
rand('state',sum(100*clock)); randn('state',sum(100*clock)); % Seeds
close all;

initPoles = 500 + 1000*(0:(nP - 1)); % Initialize pole locations
initPBW      = 80 + 40*(0:(nP - 1));
initZeros    = 600 + 1000*(0:(nZ - 1));    % Initialize with offset from pole locations by 100Hz
initZBW      = 80 + 40*(0:(nZ - 1));

initState    = [initPoles initZeros]';
initBW       = [initPBW initZBW];

% Generate data (this simply evolves pole/zeros via random walk)
[trueState BW_data] = genSynthFormantTracks(pNoiseVar, N, initState, initBW);

% Now formant tracks have been created, generate observation sequence
[y, oNoiseVar] = genNoisyObserZ(snr_dB, trueState(1:nP,:), BW_data(1:nP,:), trueState(nP+1:end,:), BW_data(nP+1:end,:), cepOrder, fs);

% Do a plot of the observations
figure;
imagesc(log(abs(y))); colorbar;
title('Cepstral Coefficients');
xlabel('Frame Number');

%%%%%%%%%%%%%%%%%%%%% Tracking Algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for tracking algorithms %%%
% x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
% y_{k}   = Hx_{k} + v_k, v_k ~ N(0, R)
% We need to set the parameters: F, Q and R
% H is obtained in the EKF via linearization about the state

F = eye(nP + nZ);   % Process Matrix F
Q = diag(var(trueState,0,2));      % Get state variance from read in tracks
R = oNoiseVar*eye(cepOrder);       % Measurement noise covariance matrix R

bwFlag = 0; % 0 - Use loaded bandwidths, 1 - Average bandwidths
bwStates = genTrackBW(bwFlag,BW_data);

% A voice activity detector is not used here in the synthetic case
formantInds = ones(N,nP + nZ);

% General Settings
algFlag = [1 1]; % Select 1 to run, 0 not to
EKF = 1; EKS = 2;

% Initial state of formant trackers
x0 = trueState(:,1);
 
countTrack = 1; % Counter for storing results
% Initialize root-mean-square error matrices:
rmse    = zeros(nP + nZ, sum(algFlag));
relRmse = zeros(nP + nZ, sum(algFlag));

% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;
    [x_estEKF x_errVarEKF] = formantTrackEKSZ(y, F, Q, R, x0, formantInds, fs, bwStates, nP, smooth);

    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'g-.'};

    % Compute and Display MSE and RMSE
    for j = 1:nP+nZ
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(N);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(N);
    end

    % Display output summary and timing information
    display(['Average EKF RMSE: ' num2str(mean(rmse(:,countTrack)))]);
    countTrack = countTrack + 1;     % Increment counter
end

% Run Extended Kalman Smoother
if algFlag(EKS)
    smooth = 1;
    [x_estEKS x_errVarEKS] = formantTrackEKSZ(y, F, Q, R, x0, formantInds, fs, bwStates, nP, smooth);

    % Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKS;
    estVar(:,:,:,countTrack) = x_errVarEKS;
    titleCell(1,countTrack+1) = {'EKS'};
    titleCell(2,countTrack+1) = {'b:'};

    %Compute and Display MSE and RMSE
    for j = 1:nP
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(N);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(N);
    end

    % Display output summary and timing information
    display(['Average EKS  RMSE: ' num2str(mean(rmse(:,countTrack)))]);
end

%Initial Plotting Variables
titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};            % Color for true state plot

% A basic plotting routine to visualize results
plotStateTracksFZ(trueState,estTracks(:,:,1),titleCell(:,[1 2]), nP);
plotStateTracksFZ(trueState,estTracks(:,:,2),titleCell(:,[1 3]), nP);