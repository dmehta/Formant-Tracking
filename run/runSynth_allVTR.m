function [] = runSynth_allVTR(snr_dB, cepOrder, fs, numFormants, trackBW, varargin)
% Test tracking algorithms using synthetic data
%
% Author: Daniel Rudoy
% Created:  02/22/2010
% Modified: 02/22/2010

% INPUT
%
% snr_dB    : How much observation noise to add
% cep_order : How many cepstal coefficients to include in the observations
% fs        : Sampling rate at which the observations are made (fake here)
% num_formants: Number of formants that we should track


addpath(genpath('../')); % Set paths
%rand('state', sum(100*clock)); randn('state',sum(100*clock)); % Set seeds
rand('state', 1); randn('state',1); % Set seeds

doRuns  = 0;
doPlots = 1;

if(doRuns)

    load '../data/VTR_Timit/allVTRTracks.mat';
    numInstances = length(DATA);              % Total number of VTR utterances

    algFlag = [0 1 1 0 0]; % Select 1 to run, 0 not to
    EKF = 1; EKS = 2; EKS_EM = 3;
    rmse    = zeros(numInstances, numFormants, sum(algFlag));
    relRmse = zeros(numInstances, sum(algFlag));
    
    % For every sentence in VTR database
    for n = 1:numIterations
        
        curVTRData = DATA{n}.vtrData';                     % Get utterance
        trueState  = 1000*curVTRData(1:numFormants, :);    % Get formant values
        BW_data    = 500*curVTRData(5:(4+numFormants),:);  % Get bandwidth values

        % Generate data according to the model

        % Now formant tracks have been created, generate observation sequence
        [y, oNoiseVar] = genNoisyObser(snr_dB, trueState, BW_data, cepOrder, fs);
        numObs = length(y); % Record number of observations

        %%%%%%%%%%%%%%%%%%%%% Tracking Algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
        % y_{k}   = H_k x_{k} + v_k, v_k ~ N(0, R)

        % We need to set the parameters: F, Q and R
        % H is obtained in the EKF via linearization about the state

        %% Process noise covariance matrix Q %%%
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

        % Initial state of formant trackers
        %initFormant = trueState(:,1);
        initFormant = 500 + 1000*(0:(numFormants - 1))';
        initBW      = 80 + 40*(0:(numFormants - 1))';

        if(trackBW)
            x0 = [initFormant', initBW']';
        else
            x0 = initFormant;
        end

        countTrack = 1; % Counter for storing results
       
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
                rmse(n,j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
                relRmse(n,j,countTrack) = (rmse(n,j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
                sDisp = sprintf('Formant %d RMSE: %2.2f ', j, rmse(n,j,countTrack));
                s = [s sDisp];
            end
            disp(s); % Show individual RMSE values
            disp('EKS Results');
            disp(sprintf('Average RMSE: %2.2f', mean(rmse(n,:,countTrack))));

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
                rmse(n,j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
                relRmse(n,j,countTrack) = (rmse(n,j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
                sDisp = sprintf('Formant %d RMSE: %2.2f ', j, rmse(n,j,countTrack));
                s = [s sDisp];
            end
            disp(s); % Show individual RMSE values
            disp('EKS_EM Results');
            disp(sprintf('Average RMSE: %2.2f', mean(rmse(n,:,countTrack))));

            countTrack = countTrack + 1;    % Increment counter
        end
        clear countTrack;
        clear estTracks;
        clear estVar;
        if(mod(n,10)==0)
            disp(sprintf('Done with example %d', n));
        end
    end
    save '../results/synthVTR_EM_EKSEM_Compare.mat';
end

if(doPlots)
    load '../results/synthVTR_EM_EKSEM_Compare.mat';
    avgDiff = mean(rmse(1:numInstances,:,1) - (rmse(1:numInstances,:,2)));
    disp(avgDiff);
    
end