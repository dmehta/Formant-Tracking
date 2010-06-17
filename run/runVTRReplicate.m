function [] = runVTRReplicate(cepOrder, fs_in, numFormants, trackBW, vtrDbNum)
% Run VTR tracker on waveforms from the VTR database
%
% Author: Daniel Rudoy
% Created:  12/13/2007
% Modified: 02/15/2010

% INPUT
% cep_order :   How many cepstal coefficients to include in the observations
% fs_in     :   Sampling rate at which the observations are made (should be 16K)
% numFormants:  Number of formants that we should track
% trackBW   :   Attempt to track bandwidths (1) or not (0)
% vtrDbNum  :   Which file number (1-516) in the VTR database?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example Usage (do NOT track bandwidths here)
% runVTRReplicate(15, 16000, 3, 0, 11);
% runVTRReplicate(15, 16000, 4, 0, 210);
% Interesting examples 200, 210
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../')); % Set paths
%addpath(genpath('../../Interspeech 2007 Code/'));
rand('state',sum(100*clock)); randn('state',sum(100*clock)); % Set seeds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Setting %%%%%%%%%%%%%%%%%%%%%%%%%

% Analysis parameters
wType = 'hamming'; % Window type
wLengthMS  = 20;   % Length of window (in milliseconds)
wOverlap   = 0.5;  % Factor of overlap of window
lpcOrder   = 12;   % Number of LPC Coefficients
peCoeff    = 0.7;  % Pre-emphasis factor

% Sampling rate
fs = 2*3500;

% Correlation control
useCorr = 1; % Flag for whether or not to use correlation
if useCorr
    diagCorr = 0; % Set to 0 if VAR(1) to 1 when 4 ind AR models
    wsTruth  = 1; % Set to 0 if to estimate matrix from VTR or to 1 for WS
end

% Voice activity detection
useVAD = 1;
if useVAD
    coastJoint = 1;    % Coast all formants jointly, or not
    quantThresh = .15; % power quantile threshold for coasting
    plotVad = 0;
end

% Select which trackers to run
algFlag = [0 1 0 0]; % Select 1 to run, 0 not to
EKF = 1; EKS = 2; EKS_EM = 3; PF = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Waveform path
dataFileName = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.wav');
% Read in audio waveform
wav = wavread(dataFileName);

% Wavesurfer file path
wsFileName   = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.FRM');
% Read in wavesurfer waveform

data = load(wsFileName);
trueStateWS = data(:,1:numFormants)';
bwDataWS    = data(:,5:1:5+numFormants-1:end)';
bwDataWS    = bwDataWS/2;

% Read in 'Truth' values from VTR database
load ../data/VTR_Timit/allVTRTracks;
curData = DATA{vtrDbNum};

% Rescale appropriately
trueStateVTR = 1000*curData.vtrData(:,1:numFormants)';    % Scaling from units
bwDataVTR    = 500*curData.vtrData(:,5:(4+numFormants))'; % that stored data was in

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
wav = wav(1:wLength*wOverlap*length(trueStateVTR));

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
y = genLPCC(wav, win, wOverlap, peCoeff, lpcOrder, cepOrder);

% Truncate at most one frame from output of LPCC call
numFrames    = size(y,2);
numObs       = numFrames;
trueStateWS  = trueStateWS(1:numFormants,1:numFrames);
bwDataWS     = bwDataWS(1:numFormants,1:numFrames);
trueStateVTR = trueStateVTR(1:numFormants,1:numFrames);
bwDataVTR    = bwDataVTR(1:numFormants,1:numFrames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Voice Activity Detection %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Tracking proceeds via:
% x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
% y_{k}   = Hx_{k} + v_k, v_k ~ N(0, R)
% To call trackers, first need to set the parameters: F, Q and R
% H is obtained within the EKF/EKS framework via linearization

%%% Set process noise covariance matrix Q %%%%
if (trackBW)
    % Get state variance from read in tracks
    Q = diag([var(trueStateWS,0,2); var(bwDataWS,0,2)]/1);
else
    % Currently is estimated from the loaded data from Wavesurfer
    Qscale = 1; % For experiments scaling this value
    % Q = diag(var(trueStateWS(:,formantInds(:,1)==1),0,2))/Qscale;
    Q = eye(numFormants)*var(trueStateWS(formantInds'==1)); % to match Interspeech 2007 code in MyTrackExp.m
end

%%% Set measurement noise covariance matrix R %%%%

% Observation noise ``variance'', should decrease as the cepstral order increases
sigExp = 1;
R      = diag(1./ones(cepOrder,1)./(((1:cepOrder).^sigExp)'));

%%% Set process evolution matrix F %%%
% Estimate cross-correlation among Formant Trajectories
if(useCorr)
    if(diagCorr)
        for dim = 1:numFormants
            ar(dim) = fitVAR(trueStateWS(dim,:)',1);
        end
        F = diag(ar);
    else
        F = fitVAR(trueStateWS(:,formantInds(:,1)==1)',1);
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
    bwFlag = 1; % 0 - Use loaded bandwidths, 1 - Average bandwidths
    bwStates = genTrackBW(bwFlag, bwDataWS);
end

%%% Set initial state of formant trackers
%initFormant = trueStateWS(:,1); % Uncomment to use the first True state
% Using values suggested by Li Deng and otherwise in literature
initFormant = 500 + 1000*(0:(numFormants - 1))';
initBW      = 80 + 40*(0:(numFormants - 1))';
if(trackBW)
    x0 = [initFormant', initBW']';
else
    x0 = initFormant;
end
display(['Initial state estimates set at: ' num2str(x0')])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countTrack = 1; % Counter for storing results
% Initialize root-mean-square error matrices:
rmse    = zeros(numFormants, sum(algFlag));
relRmse = zeros(numFormants, sum(algFlag));

% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;  % No smoothing
    [x_estEKF x_errVarEKF] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth);

    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'g-.'};

    % Compute and Display MSE and RMSE
    for j = 1:numFormants
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueStateVTR(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR(j,:)))*sqrt(numObs);
    end
    display(['Average RMSE: ' num2str(mean(rmse(:,countTrack)))]);
    countTrack = countTrack + 1;     % Increment counter
end

% Run Extended Kalman Smoother
if algFlag(EKS)
    smooth = 1;
    %Q = 1.0458e+006*eye(numFormants);
    %Q = 8.0244e+005*eye(numFormants);
    [x_estEKS x_errVarEKS] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth);
      
    % Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKS;
    estVar(:,:,:,countTrack) = x_errVarEKS;
    titleCell(1,countTrack+1) = {'EKS'};
    titleCell(2,countTrack+1) = {'b:'};

    %Compute and Display MSE and RMSE
    s = [];
    for j = 1:numFormants
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueStateVTR(j,:)))/sqrt(numObs);

        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR(j,:)))*sqrt(numObs);

        % Compute errors where there was speech energy (i.e., omit silences)
        Sinds = find(formantInds(:,j) == 1);
        rmseS(j,countTrack) = norm((estTracks(j,Sinds,countTrack)-trueStateVTR(j,Sinds)))/sqrt(numObs);
        relRmseS(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR(j,Sinds)))*sqrt(numObs);

        sDisp  = sprintf('Formant %d RMSE: %2.2f; NS RMSE: %2.2f', j, rmse(j,countTrack), rmseS(j,countTrack));
        s = [s sDisp];
    end
    disp('Extended Kalman Smoother Results');
    disp(s); % Show individual RMSE values
    disp(sprintf('Average RMSE: %2.2f; NS RMSE: %2.2f', mean(rmse(:,countTrack)),mean(rmseS(:,countTrack))));
    countTrack = countTrack + 1;
end

% Run Extended Kalman Smooother with EM steps
if algFlag(EKS_EM)

    maxNumIter = 1;  % Maximum number of EKS-EM iterations

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

    % Compute and display RMSE and relative RMSE
    s = [];
    for j = 1:numFormants
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueStateVTR(j,:)))/sqrt(numObs);

        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR(j,:)))*sqrt(numObs);

        % Compute errors where there was speech energy (i.e., omit silences)
        Sinds = find(formantInds(:,j) == 1);
        rmseS(j,countTrack) = norm((estTracks(j,Sinds,countTrack)-trueStateVTR(j,Sinds)))/sqrt(numObs);
        relRmseS(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR(j,Sinds)))*sqrt(numObs);

        sDisp  = sprintf('Formant %d RMSE: %2.2f; NS RMSE: %2.2f', j, rmse(j,countTrack), rmseS(j,countTrack));
        s = [s sDisp];
    end
    disp('Extended Kalman Smoother Results');
    disp(s); % Show individual RMSE values
    disp(sprintf('Average RMSE: %2.2f; NS RMSE: %2.2f', mean(rmse(:,countTrack)),mean(rmseS(:,countTrack))));
    countTrack = countTrack + 1;    % Increment counter
end


if algFlag(PF)
    numParticles = 1000;
    [x_estEKS x_errVarEKS] = formantTrackPF(y, F, Q, R, x0, formantInds, fs, bwStates, numParticles);

    % Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKS;
    %estVar(:,:,:,countTrack) = x_errVarEKS;
    titleCell(1,countTrack+1) = {'PF'};
    titleCell(2,countTrack+1) = {'b:'};

    %Compute and Display MSE and RMSE
    s = [];
    for j = 1:numFormants
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueStateVTR(j,:)))/sqrt(numObs);

        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR(j,:)))*sqrt(numObs);

        % Compute errors where there was speech energy (i.e., omit silences)
        Sinds = find(formantInds(:,j) == 1);
        rmseS(j,countTrack) = norm((estTracks(j,Sinds,countTrack)-trueStateVTR(j,Sinds)))/sqrt(numObs);
        relRmseS(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR(j,Sinds)))*sqrt(numObs);

        sDisp  = sprintf('Formant %d RMSE: %2.2f; NS RMSE: %2.2f', j, rmse(j,countTrack), rmseS(j,countTrack));
        s = [s sDisp];
    end
    disp('Particle Filter Results');
    disp(s); % Show individual RMSE values
    disp(sprintf('Average RMSE: %2.2f; NS RMSE: %2.2f', mean(rmse(:,countTrack)),mean(rmseS(:,countTrack))));
    countTrack = countTrack + 1;
end


%Initial Plotting Variables
titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};            % Color for true state plot

% A basic plotting routine to visualize results
plotStateTracks(trueStateVTR,estTracks,titleCell);

% Super-impose over a spectrogram

plotSpecTracks2(wav, estTracks(:,:,1), aParams);
gcf;
title('EKS');

plotSpecTracks2(wav, trueStateVTR, aParams);
gcf;
title('VTR');

