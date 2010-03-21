function [] = runSynthZ(testMethod, snr_dB, cepOrder, fs, numFormants, numZeros, varargin)

% Track poles and zeros (no bandwidths) and on synthetic data
% Author: Daniel Rudoy
% Created:  12/13/2007
% Modified: 12/13/2007, 02/15/2010

% INPUT
%
% testMethod: Synth, VTR, Praat, WS
% snr_dB    : How much observation noise to add
% cep_order : How many cepstal coefficients to include in the observations
% fs        : Sampling rate at which the observations are made (fake here)
% num_formants: Number of formants that we should track

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% Synthetic Example:
% runSynthZ('Synth',15, 15, 16000, 4, 4, 100, 30^2);
%
% The usual examples from PRAAT and Wavesurfer are not included, because
% it is not clear how to generate paths of zeros along with the paths of
% formants that are available.

addpath(genpath('../')); % Paths
rand('state',sum(100*clock)); randn('state',sum(100*clock)); % Seeds
close all;

if(length(varargin) ~= 2)
    error('Please provide number formants, number observations and the variance of the process noise');
else
    numObs      = varargin{1}; % Number of observations to generate
    pNoiseVar   = varargin{2}; % Variance of the process noise
end


initFormants = 500 + 1000*(0:(numFormants - 1));
initFBW      = 80 + 40*(0:(numFormants - 1));
initZeros    = 550 + 1000*(0:(numZeros - 1)); % Initialize with offset from formants
initZBW      = 80 + 40*(0:(numZeros - 1));

initState    = [initFormants initZeros]';
initBW       = [initFBW initZBW];

% Generate data (this simply evolves pole/zeros via random walk)
[trueState BW_data] = genSynthFormantTracks(pNoiseVar, numObs, initState, initBW);

% Now formant tracks have been created, generate observation sequence
[y, oNoiseVar] = genNoisyObserZ(snr_dB, trueState(1:numFormants,:), BW_data(1:numFormants,:), trueState(numFormants+1:end,:), BW_data(numFormants+1:end,:), cepOrder, fs);

numObs = length(y); % Record number of observations

% Do a plot of the observations
figure;
imagesc(log(abs(y)));
colorbar;
title('Cepstral Coefficients');

%%%%%%%%%%%%%%%%%%%%% Tracking Algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for tracking algorithms %%%
% x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
% y_{k}   = Hx_{k} + v_k, v_k ~ N(0, R)
% We need to set the parameters: F, Q and R
% H is obtained in the EKF via linearization about the state

F = eye(numFormants + numZeros);   % Process Matrix F
Q = diag(var(trueState,0,2));      % Get state variance from read in tracks
R = oNoiseVar*eye(cepOrder);       % Measurement noise covariance matrix R

bwFlag = 0; % 0 - Use loaded bandwidths, 1 - Average bandwidths
bwStates = genTrackBW(bwFlag,BW_data);

% A voice activity detector is not used here in the synthetic case
formantInds = ones(numObs,numFormants+numZeros);

% General Settings
algFlag = [1 1]; % Select 1 to run, 0 not to
EKF = 1; EKS = 2;

% Initial state of formant trackers
x0 = trueState(:,1);
 
countTrack = 1; % Counter for storing results
% Initialize root-mean-square error matrices:
rmse    = zeros(numFormants + numZeros, sum(algFlag));
relRmse = zeros(numFormants + numZeros, sum(algFlag));

% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;
    [x_estEKF x_errVarEKF] = formantTrackEKSZ(y, F, Q, R, x0, formantInds, fs, bwStates, numFormants,0);

    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'g-.'};

    % Compute and Display MSE and RMSE
    for j = 1:numFormants+numZeros
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
    end

    % Display output summary and timing information
    display(['Average RMSE: ' num2str(mean(rmse(:,countTrack)))]);
    countTrack = countTrack + 1;     % Increment counter
end

% Run Extended Kalman Smoother
if algFlag(EKS)
    smooth = 1;
    [x_estEKS x_errVarEKS] = formantTrackEKSZ(y, F, Q, R, x0, formantInds, fs, bwStates, numFormants, smooth);

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
    display(['Average RMSE: ' num2str(mean(rmse(:,countTrack)))]);
end

%Initial Plotting Variables
titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};            % Color for true state plot

% A basic plotting routine to visualize results
plotStateTracks(trueState,estTracks(:,:,1),titleCell(:,[1 2]));
plotStateTracks(trueState,estTracks(:,:,2),titleCell(:,[1 3]));
