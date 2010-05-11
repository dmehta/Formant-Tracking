function [rmse] = runSynth(testMethod, snr_dB, cepOrder, fs, numFormants, trackBW, numParticles, doPlots, varargin)
% Test tracking algorithms using synthetic data
%
% Author: Daniel Rudoy
% Created:  12/13/2007
% Modified: 02/15/2010

% INPUT
%
% testMethod    : Synth, VTR, PRAAT, WS
% snr_dB        : How much observation noise to add
% cepOrder     : How many cepstal coefficients to include in the observations
% fs            : Sampling rate at which the observations are made (fake here)
% numFormants   : Number of formants that we should track
% trackBW       : Flag whether to track (1) or not track (0) bandwidths
% numParticles  : ?
% doPlots       : Flag for whether to plot (1) or not to plot (0) estimates
%                       and ground truth tracks
% varargin      : Dependent on testMethod. If
%                       Synth: number of observations, process noise variance
%                       VTR: data filename, sample index
%                       PRAAT, WS: data filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% Synthetic Example:
% runSynth('Synth',105, 15, 16000, 4, 0, [], 1, 30, 2);
%
% VTR Database Example:
% runSynth('VTR', 15, 15, 16000, 4, 1, [], 1, '../data/VTR_Timit/allVTRTracks.mat', 35);
%
% PRAAT and Wavesurfer Examples. PRAAT has .Formant output files, while
% Wavesurfer has .FRM output files:
% runSynth('PRAAT',105, 15, 10000, 4, [], 0, 1, '../data/synthData/ln.roy.10k.Formant');
% runSynth('PRAAT',105, 15, 10000, 4, [], 0, 1, '../data/synthData/mlm.tea.10k.Formant');
% runSynth('WS', 15, 15, 10000, 4, 1, [], 1, '../data/synthData/ln.roy.10k.FRM');
% runSynth('WS', 20, 15, 10000, 4, 0, [], 1, '../data/synthData/mlm.tea.10k.FRM');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all;

% Error checking on input
switch testMethod

    case 'Synth'
        if(length(varargin) ~= 2)
            error('Please provide number observations and the variance of the process noise');
        else
            numObs      = varargin{1}; % Number of observations to generate
            pNoiseVar   = varargin{2}; % Variance of the process noise
        end

    case {'VTR'}
        if(length(varargin) ~= 2)
            error('Check parameter usage: data filename and sample index');
        else
            dataFileName = varargin{1};
            sampleIndex  = varargin{2};
        end

    case{'PRAAT', 'WS'}
        if(length(varargin) ~= 1)
            error('Must provide filename from which to read formant tracks');
        else
            dataFileName = varargin{1};
        end

    otherwise
        error('Invalid Test Command. Allowed: Synth, VTR, WS, PRAAT');
end

addpath(genpath('../')); % Set paths
rand('state',sum(100*clock)); randn('state',sum(100*clock)); % Seeds
% rand('state', 1); randn('state',1); % Set seeds
% doPlots = 1;

% Generate data according to the model
switch testMethod
    case 'Synth' % Synthetic Data Test
        display('Generating formant data');
        % Number of observations
        initState = 500 + 1000*(0:(numFormants - 1))'; % Note: matched to track init
        %initState = [500 1500 2500 3500];
        initBW    = 80 + 40*(0:(numFormants - 1));     % Note: matched to track init

        % Generate Data
        [trueState BW_data] = genSynthFormantTracks(pNoiseVar, numObs, initState, initBW);
    case {'VTR', 'PRAAT', 'WS'} % VTR or other comparative data test
        switch testMethod
            case 'VTR'
                display(['Applying to example ' int2str(sampleIndex) ' of VTR database file ' dataFileName])
                [trueState, BW_data] = vtrFormantRead(dataFileName, numFormants, sampleIndex);
            case 'PRAAT'
                display(['Reading data from Praat-generated file ' dataFileName])
                [trueState, BW_data] = praatFormantRead(dataFileName,numFormants);
            case 'WS'
                display(['Reading data from Wavesurfer-generated file ' dataFileName])
                [trueState, BW_data] = wavesurferFormantRead(dataFileName,numFormants);
        end
end

% Now that formant tracks have been created, generate observation sequence
[y, oNoiseVar] = genNoisyObser(snr_dB, trueState, BW_data, cepOrder, fs);

numObs = length(y); % Record number of observations

rand('state', sum(100*clock)); randn('state',sum(100*clock)); % Set seeds


%%%%%%%%%%%%%%%%%%%%% Tracking Algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
% y_{k}   = Hx_{k} + v_k, v_k ~ N(0, R)

% We need to set the parameters: F, Q and R
% H is obtained in the EKF via linearization about the state

%%% Process noise covariance matrix Q %%%%
if(trackBW)
    %Q = diag([var(trueState,0,2); var(BW_data,0,2)]);            % Get state variance from read in tracks
    Q = diag([var(trueState(:,2:end)-trueState(:,1:end-1),0,2); var(BW_data(:,2:end) - BW_data(:,1:end-1),0,2); ]);
    F = eye(2*numFormants);       % Process matrix F
else
    %Q = diag(var(trueState,0,2)); % Get state variance from read in tracks
    Q = diag(var(trueState(:,2:end)-trueState(:,1:end-1),0,2));
    F = eye(numFormants);         % Process matrix F
end

R = oNoiseVar*eye(cepOrder);  % Measurement noise covariance matrix R

%%% Set Bandwidths %%%%
if(trackBW)
    bwStates = []; % If we are tracking bandwidths do not provide them
else
    bwFlag = 0; % 0 - Use loaded bandwidths, 1 - Average bandwidths
    bwStates = genTrackBW(bwFlag,BW_data);
end

% A voice activity detector is not used here in the synthetic case
formantInds = ones(numObs,numFormants);

% General Settings
% Decision Flags for Tracking Processes, and parameters
algFlag = [1 0 0 0 0]; % Select 1 to run, 0 not to
EKF = 1; EKS = 2; EKS_EM = 3; PF = 4; RBPF = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial state of formant trackers
%initFormant = trueState(:,1);
initFormant = 500 + 1000*(0:(numFormants - 1))';
initBW = 80 + 40*(0:(numFormants - 1))';

if(trackBW)
    x0 = [initFormant', initBW']';
else
    x0 = initFormant;
end

display(['Initial state estimates set at: ' num2str(x0')])

countTrack = 1; % Counter for storing results

% Initialize root-mean-square error matrices:
if(trackBW)
    rmse    = zeros(numFormants*2, sum(algFlag));
    relRmse = zeros(numFormants*2, sum(algFlag));
else
    rmse    = zeros(numFormants, sum(algFlag));
    relRmse = zeros(numFormants, sum(algFlag));
end

% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;
    [x_estEKF x_errVarEKF] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth);

    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'g-.'};

    % Compute and display RMSE and relative RMSE
    if(trackBW)
        trueStateF_BW = [trueState; BW_data];
        for j = 1:numFormants*2
            rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueStateF_BW(j,:)))/sqrt(numObs);
            relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueStateF_BW(j,:)))*sqrt(numObs);
        end
        
    else
        for j = 1:numFormants
            rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
            relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
        end
    end        

    display(['Average RMSE: ' num2str(mean(rmse(:,countTrack)))]);
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
    titleCell(2,countTrack+1) = {'b-'};

    % Compute and display RMSE and relative RMSE
    s = [];
    for j = 1:numFormants
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
        sDisp = sprintf('Formant %d RMSE: %2.2f ', j, rmse(j,countTrack));
        s = [s sDisp];
    end
    disp(s); % Show individual RMSE values
    disp('EKS Results');
    disp(sprintf('Average RMSE: %2.2f', mean(rmse(:,countTrack))));

    countTrack = countTrack + 1;    % Increment counter
end

if algFlag(EKS_EM)

    % Maximum number of EKS-EM iterations
    maxNumIter = 10;
    % Load initial values
    thetaInit.F = F;
    thetaInit.Q = Q;
    thetaInit.R = R;
    thetaInit.x0 = x0;

    [m_upS P_upS theta logL] = formantTrackEKS_EM(y, thetaInit, formantInds, fs, bwStates, maxNumIter);

    % Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = m_upS;
    estVar(:,:,:,countTrack)  = P_upS;
    titleCell(1,countTrack+1) = {'EKS EM'};
    titleCell(2,countTrack+1) = {'k-'};

    s = [];
    for j = 1:numFormants
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
        sDisp = sprintf('Formant %d RMSE: %2.2f ', j, rmse(j,countTrack));
        s = [s sDisp];
    end
    disp(s); % Show individual RMSE values
    disp('EKS_EM Results');
    disp(sprintf('Average RMSE: %2.2f', mean(rmse(:,countTrack))));

    countTrack = countTrack + 1;    % Increment counter
end

% Run Particle Filter
if algFlag(PF)

    t0 = clock;
    [pmean phist Neff] = formantTrackPF(y, F, Q, R, x0, formantInds, fs, BW_data, numParticles);

    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = pmean;
    titleCell(1,countTrack+1) = {'PF'};
    titleCell(2,countTrack+1) = {'c--'};

    %Compute and Display MSE and RMSE
    for j = 1:numFormants
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
    end

    % Display output summary and timing information
    display([int2str(numParticles) '-Particle Filter Run Time: ' num2str(etime(clock,t0)) ' s. Average RMSE: ' num2str(mean(rmse(:,countTrack)))]);

    % Increment counter
    countTrack = countTrack + 1;
end

% Execute Rao-Blackwellized Particle Filter
if algFlag(RBPF)

    t0 = clock;
    formantInds = ones(size(y,2),numFormants);
    numParticles = 20;
    [x_estRBPF,x_errVarRBPF,pEst] = formantTrackRBPF(numParticles, y, formantInds, pNoiseVar, oNoiseVar,fs,trBW_flag,BW_data,initState,F);

    disp(strcat('MMSE estimate of lambda: ', num2str(1/pEst)));
    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estRBPF;
    %estVar(:,:,:,countTrack) = x_errVarRBPF; FIX THIS LATER
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'k-.'};

    %Compute and Display MSE and RMSE
    for j = 1:numFormants
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
    end

    % Display output summary and timing information
    display(['RBPF Run Time: ' num2str(etime(clock,t0)) ' s. Average RMSE: ' num2str(mean(rmse(:,countTrack)))]);
end

if(doPlots)

    %Initial Plotting Variables
    titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
    titleCell(2,1)  = {'r'};            % Color for true state plot

    % A basic plotting routine to visualize results
    if(trackBW)
        plotStateTracks([trueState; BW_data],estTracks,titleCell);
    else
        for j = 1:sum(algFlag)
            plotStateTracks(trueState,estTracks(:,:,j),titleCell(:,[1 j+1]));
        end
    end
end