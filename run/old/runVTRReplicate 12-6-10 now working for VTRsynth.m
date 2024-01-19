% addpath(genpath('../')); % Set paths

%%
cepOrder = 20;
fs = 2*8000;
useCorr = 0; % Flag for whether or not to use correlation
useVAD = 0; % 0.3 s, all disp 0.6 s
numFormants = 4; % to track
trackBW = 0;
vtrDbNum = 265;
zOrder = 0;

% Select which trackers to run
algFlag = [0 1 0 0]; % Select 1 to run, 0 not to
EKF = 1; EKS = 2; EKS_EM = 3; PF = 4;

% function [] = runVTRReplicate(cepOrder, fs_in, numFormants, trackBW, vtrDbNum)
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

rand('state',sum(100*clock)); randn('state',sum(100*clock)); % Set seeds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Setting %%%%%%%%%%%%%%%%%%%%%%%%%

% Analysis parameters
wType = 'hamming'; % Window type
wLengthMS  = 20;   % Length of window (in milliseconds)
wOverlap   = 0.5;  % Factor of overlap of window
lpcOrder   = 12;   % Number of LPC Coefficients
peCoeff    = 0.7;  % Pre-emphasis factor

% Correlation control
if useCorr
    diagCorr = 0; % Set to 0 if VAR(1) to 1 when 4 ind AR models
    wsTruth  = 1; % Set to 0 if to estimate matrix from VTR or to 1 for WS
end

% Voice activity detection
if useVAD
    coastJoint = 1;    % Coast all formants jointly, or not
    quantThresh = .15; % power quantile threshold for coasting
    plotVad = 0;
end

rmseAll = zeros(516, 3, 4); % database length X num algorithms (EKS, WS, Praat) X formant number (1, 2, 3, overall)
rmseAllS = zeros(516, 3, 4); % same dims, but speech frames only

%%
% tic
for ii = 1%:516 % praat bombs on 364
    %%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vtrDbNum = ii;
    %disp('Database number')
    disp(vtrDbNum)
    
    % Waveform path
    dataFileName = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.wav');
    dataFileNameIn = strcat('../data/VTRsynth/VTRsynth',num2str(vtrDbNum),'.wav');
    
    % Read in audio waveform
    [wav_orig, fs_in] = wavread(dataFileNameIn);

    % Wavesurfer file path
    wsFileName   = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.FRM');
    % Read in wavesurfer waveform

    data = load(wsFileName);
    trueStateWS = data(:,1:numFormants)';
    bwDataWS    = data(:,numFormants+1:2*numFormants)';
    bwDataWS    = bwDataWS/2; % why divide by 2? *****************

    % Read in 'Truth' values from VTR database
    load ../data/VTR_Timit/allVTRTracks;
    curData = DATA{vtrDbNum};
    clear DATA

    % Rescale appropriately
    trueStateVTR = 1000*curData.vtrData(:,1:numFormants)';    % Scaling from units
    bwDataVTR    = 500*curData.vtrData(:,5:(4+numFormants))'; % that stored data was in, why not *1000? ******************

    %%%%%%%%%%%%%%%%%%%%%% Generate observation sequence %%%%%%%%%%%%%%%%%%%%%%
    % Resample input data if sampling rate is not equal to that of input
    if(fs ~= fs_in)
        wav = resample(wav_orig,fs,fs_in,2048);
        %display(['Input Fs = ' int2str(fs_in) ' Hz; resampling to ' num2str(fs) ' Hz'])
    else
        wav = wav_orig;
    end

%     trueStateVTR = trueStateVTR(:, 1:numFrames);
%     bwDataVTR = bwDataVTR(:, 1:numFrames);
%     trueStateWS = trueStateWS(:, 1:numFrames);
%     bwDataWS = bwDataWS(:, 1:numFrames);
    
    % Compute window length in samples, now that sampling rate is set
    wLength = floor(wLengthMS/1000*fs);
    wLength = wLength + (mod(wLength,2)~=0); % Force even
    win = feval(wType,wLength);

    % Truncate the waveform to match true state
%     wav = wav(1:wLength*wOverlap*(length(trueStateVTR)+1));

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
    % y = genLPCC(wav, win, wOverlap, peCoeff, lpcOrder, cepOrder);
    y = genLPCCz(wav, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);

    % Truncate at most one frame from output of LPCC call
    numFrames    = min([size(y,2), size(trueStateWS, 2), size(trueStateVTR, 2)]);
    numObs       = numFrames;
    
    y            = y(:, 1:numFrames);
    trueStateWS  = trueStateWS(1:numFormants,1:numFrames);
    bwDataWS     = bwDataWS(1:numFormants,1:numFrames);
    trueStateVTR = trueStateVTR(1:numFormants,1:numFrames);
    bwDataVTR    = bwDataVTR(1:numFormants,1:numFrames);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Speech Activity Detection %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(useVAD)
        %display('Using Voice Activity Detection');
        % If all formants are coasted jointly, then do not do multiband detection
        multiBand = ~coastJoint;
        
        % simple VAD
        %[frameInds] = simpleVAD(wav, aParams, quantThresh, multiBand, [], plotVad);
        
        % TIMIT labels to help define speech frames
        load([dataFileName(1:end-4), '.mat']) % 0.44 s
        frameInds = timitVAD(wav, data.phones, aParams, plotVad); % 0.3 s
        
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
            if trackBW
                formantInds = repmat(frameInds,2*numFormants,1)';
            else
                formantInds = repmat(frameInds,numFormants,1)';
            end
        end
    else
        %display('Not using Voice Activity Detection');
        if trackBW
            formantInds = ones(numFrames,2*numFormants);
        else
            formantInds = ones(numFrames,numFormants);
        end
    end

    % get rid of indices that are out of range
    formantInds = formantInds(1:numFrames, :);
    
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
    %R = 1/100*diag(1./ones(cepOrder,1));

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
    %display(['Initial state estimates set at: ' num2str(x0')])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    countTrack = 1; % Counter for storing results
    % Initialize root-mean-square error matrices:
    rmse    = zeros(numFormants, sum(algFlag));
    relRmse = zeros(numFormants, sum(algFlag));
    clear estTracks estVar

    % Run Extended Kalman Filter
    if algFlag(EKF)
        smooth = 0;  % No smoothing
        %[x_estEKF x_errVarEKF] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth);
        [x_estEKF x_errVarEKF] = formantTrackEKSZ(y, F, Q, R, x0, formantInds', fs, bwStates, numFormants, smooth);

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
        %display(['Average RMSE: ' num2str(mean(rmse(:,countTrack)))]);
        countTrack = countTrack + 1;     % Increment counter
    end

    % Run Extended Kalman Smoother
    if algFlag(EKS)
        smooth = 1;
        %Q = 1.0458e+006*eye(numFormants);
        %Q = 8.0244e+005*eye(numFormants);
        %[x_estEKS x_errVarEKS] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth);
        [x_estEKS x_errVarEKS] = formantTrackEKSZ(y, F, Q, R, x0, formantInds', fs, bwStates, numFormants, smooth);

        % Track estimate into data cube for plot routines
        estTracks(:,:,countTrack) = x_estEKS;
        estVar(:,:,:,countTrack) = x_errVarEKS;
        titleCell(1,countTrack+1) = {'EKS'};
        titleCell(2,countTrack+1) = {'b:'};

%         %Compute and Display MSE and RMSE
%         s = [];
%         for j = 1:numFormants
%             rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueStateVTR(j,:)))/sqrt(numObs);
% 
%             relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR(j,:)))*sqrt(numObs);
% 
%             % Compute errors where there was speech energy (i.e., omit silences)
%             Sinds = find(formantInds(:,j) == 1);
%             rmseS(j,countTrack) = norm((estTracks(j,Sinds,countTrack)-trueStateVTR(j,Sinds)))/sqrt(numObs);
%             relRmseS(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR(j,Sinds)))*sqrt(numObs);
% 
%             sDisp  = sprintf('Formant %d RMSE: %2.2f; NS RMSE: %2.2f', j, rmse(j,countTrack), rmseS(j,countTrack));
%             s = [s sDisp];
%         end
%         disp('Extended Kalman Smoother Results');
%         disp(s); % Show individual RMSE values
%         disp(sprintf('Average RMSE: %2.2f; NS RMSE: %2.2f', mean(rmse(:,countTrack)),mean(rmseS(:,countTrack))));
%         countTrack = countTrack + 1;
    end

    % Run Extended Kalman Smooother with EM steps
    if algFlag(EKS_EM)

        maxNumIter = 100;  % Maximum number of EKS-EM iterations

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
    
    %% find RMSE
    countTrack = 1;
    rmse    = zeros(numFormants, sum(algFlag));
    relRmse = zeros(numFormants, sum(algFlag));
    rmseS    = zeros(numFormants, sum(algFlag));
    relRmseS = zeros(numFormants, sum(algFlag));    
    
    % choose tracker: 1. EKS, 2. WaveSurfer, 3. Praat
    for jj = 1:3
        trackerType = jj;

        switch trackerType
            case 1  
                %len = min(size(trueStateVTR,2),size(estTracks,2));
                %estTracks2 = estTracks(:, 1:len, countTrack);
                estTracks2 = estTracks;
                trueStateVTR2 = trueStateVTR;
                algo = 'Extended Kalman Smoother Results';
            case 2
                %[f,bw] = wavesurferFormantRead([wsFileName(1:end-3) 'frm'],numFormants);
                %f = trueStateWS;
                %len = min(size(trueStateVTR,2),size(f,2));
                %estTracks2 = f(:, 1:len);
                estTracks2 = trueStateWS;
                trueStateVTR2 = trueStateVTR;
                algo = 'WaveSurfer Results';
            case 3
                praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
                params = getFormantParamsDefault;
                 params.numTracks = numFormants;
                 params.winLen = aParams.wLengthMS/1000;
                 params.timestep = params.winLen*aParams.wOverlap;
                 params.maxformantHz = fs/2;
                 params.peHz = params.maxformantHz;
                 params.maxformantcount = numFormants + 2; % for some reason this has to be done for Praat to work
                [f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(wav, aParams.fs, praat_dir, params);
                f = f';
                % Praat's first values at 20 ms, not 10 ms like VTR/TIMIT
                %len = min(size(trueStateVTR,2),size(f,2));
                %estTracks2 = f(:, 1:len);
                trueStateVTR2 = trueStateVTR(:, 2:end, countTrack);
                len = min(size(trueStateVTR2, 2), size(f, 2));
                
                trueStateVTR2 = trueStateVTR2(:, 1:len);
                estTracks2 = f(:, 1:len);
                formantInds = formantInds(1:len, :);
%                 if size(estTracks2, 2) < size(trueStateVTR2, 2)
%                     trueStateVTR2 = trueStateVTR2(:, 1:end-1, countTrack); % edge effect
%                 end
                algo = 'Praat Results';
            otherwise
                error('Unknown tracker type.')
        end

        %Compute and Display MSE and RMSE
        %trueStateVTR2 = trueStateVTR;
        s = [];
        for j = 1:numFormants
            rmse(j,countTrack) = norm((estTracks2(j,:,countTrack)-trueStateVTR2(j,:)))/sqrt(numObs);
            relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR2(j,:)))*sqrt(numObs);

            % Compute errors where there was speech energy (i.e., omit silences)
            Sinds = find(formantInds(:,j) == 1);
            rmseS(j,countTrack) = norm((estTracks2(j,Sinds,countTrack)-trueStateVTR2(j,Sinds)))/sqrt(numObs);
            relRmseS(j,countTrack) = (rmse(j,countTrack)/norm(trueStateVTR2(j,Sinds)))*sqrt(numObs);

            rmseAll(vtrDbNum, trackerType, j) = rmse(j,countTrack);
            rmseAllS(vtrDbNum, trackerType, j) = rmseS(j,countTrack);
            
            sDisp  = sprintf('Formant %d RMSE: %2.2f; NS RMSE: %2.2f\n', j, rmse(j,countTrack), rmseS(j,countTrack));
            s = [s sDisp];            
        end
        rmseAll(vtrDbNum, trackerType, numFormants+1) = mean(rmse(:,countTrack));
        rmseAllS(vtrDbNum, trackerType, numFormants+1) = mean(rmseS(:,countTrack));
        
        % Show individual RMSE values
        disp(algo);
        fprintf(s);
        disp(sprintf('Average   RMSE: %2.2f; NS RMSE: %2.2f\n', mean(rmse(:,countTrack)),mean(rmseS(:,countTrack))));        
    end
    disp(' ')
end
% toc
% save('..\results\rmse_EKS_WS_Praat_trackBW0_cep12_ar12ma0_win20_pe07_fs7000_simpleVAD_useCorr1', 'rmseAll', 'rmseAllS')
% save('..\results\rmse_EKS_WS_Praat_trackBW0_cep12_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1', 'rmseAll', 'rmseAllS')
% save('..\results\rmse_EKS_WS_Praat_trackBW1_cep12_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1', 'rmseAll', 'rmseAllS')
% save('..\results\rmse_EKS_WS_Praat_trackBW0_cep12_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr0', 'rmseAll', 'rmseAllS')
% save('..\results\rmse_EKS_WS_Praat_trackBW1_cep12_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr0')

% save('..\results\rmseVTRsynth_EKS_WS_Praat_trackBW0_cep12_ar12ma0_win20_pe07_fs16000_timitVAD0_useCorr0_db1-363')
% save('..\results\rmseVTRsynth_EKS_WS_Praat_trackBW1_cep12_ar12ma0_win20_pe07_fs16000_timitVAD0_useCorr0_db1-363')

break

%%
% load('..\results\rmse_EKS_WS_Praat_trackBW0_cep12_ar12ma0_win20_pe07_fs7000_simpleVAD_useCorr1', 'rmseAll', 'rmseAllS')
% load('..\results\rmse_EKS_WS_Praat_trackBW0_cep12_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1', 'rmseAll', 'rmseAllS')
% load('..\results\rmse_EKS_WS_Praat_trackBW1_cep12_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1', 'rmseAll', 'rmseAllS')
% load('..\results\rmse_EKS_WS_Praat_trackBW0_cep12_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr0', 'rmseAll', 'rmseAllS')
% load('..\results\rmse_EKS_WS_Praat_trackBW1_cep12_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr0', 'rmseAll', 'rmseAllS')
% load('..\results\rmseVTRsynth_EKS_WS_Praat_trackBW0_cep12_ar12ma0_win20_pe07_fs16000_timitVAD0_useCorr0_db1-363')
% load('..\results\rmseVTRsynth_EKS_WS_Praat_trackBW1_cep12_ar12ma0_win20_pe07_fs16000_timitVAD0_useCorr0_db1-363')

for ii = 1:3
    for jj = 1:4
        %tmp = rmseAllS([1:264 266:end], ii, jj);
        tmp = rmseAllS(:, ii, jj);
        %rmse_table(ii, jj) = mean(tmp(~isnan(tmp))); % take into account any NaNs
        rmse_table(ii, jj) = mean(tmp(1:363));
    end
end

figure, plot(rmseAllS(1:end, 1, 1))

rmse_table

%% A basic plotting routine to visualize results
plotStateTracks(trueStateVTR,estTracks,titleCell);

%% additional plotting routines Nov 2010
figure, plot(trueStateVTR')
hold on, plot(bwDataVTR')
% figure, plot(trueStateWS')
figure, plot(estTracks')

%% state tracks with uncertainties and truth
plotStateTracksFZ_EstVar_Truth(trueStateVTR,estTracks,estVar,titleCell,numFormants,trackBW)
disp(['Mean RMSE: ', num2str(mean(rmse))])
rmse

%% Super-impose Ground Truth over a spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
numAntiF = 0;

plotSpecTracks2BW(wav, trueStateVTR, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Super-impose EKS over a spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
numAntiF = 0;

plotSpecTracks2BW(wav, estTracks, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Wavesurfer tracks on spectrogram
% [f,bw] = wavesurferFormantRead([wsFileName(1:end-3) 'frm'],numFormants);
% plotSpecTracksWS(wav, f, bw, aParams, trackBW)
plotSpecTracksWS(wav, trueStateWS, bwDataWS, aParams, trackBW)
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
set(gca, 'PlotBoxAspectRatio', [10 1 1])
 
%% Praat tracks on spectrogram
praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
params = getFormantParamsDefault;
 params.winLen = aParams.wLengthMS/1000;
 params.timestep = params.winLen*aParams.wOverlap;
 params.numTracks = numFormants;
 params.maxFormantHz = fs/2;
 params.peHz = params.maxformantHz;
[f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(wav, aParams.fs, praat_dir, params);
plotSpecTracksPraat(wav, f', bw', aParams, formantTimes, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and EKS on spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
numAntiF = 0;

plotSpecTracks2BW_Truth(wav, trueStateVTR, estTracks, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and WaveSurfer on spectrogram
[f,bw] = wavesurferFormantRead([wsFileName(1:end-3) 'frm'],numFormants);
plotSpecTracksWS_Truth(wav, f, bw, trueStateVTR, aParams, trackBW)
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and Praat tracks on spectrogram
praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
params = getFormantParamsDefault;
 params.numTracks = numFormants;
 params.winLen = aParams.wLengthMS/1000;
 params.timestep = params.winLen*aParams.wOverlap;
 params.maxformantHz = fs/2;
 params.peHz = params.maxformantHz;
 params.maxformantcount = numFormants + 2; % for some reason this has to be done for Praat to work
[f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(wav, aParams.fs, praat_dir, params);
plotSpecTracksPraat_Truth(wav, f', bw', trueStateVTR, aParams, formantTimes, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% EKS with uncertainty on spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
numAntiF = 0;

plotSpecTracks2_estVar(wav, estTracks, x_errVarEKS, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and EKS with uncertainty on spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
numAntiF = 0;

plotSpecTracks2_estVar_Truth(wav, estTracks, x_errVarEKS, trueStateVTR, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and EKS with uncertainty on spectrogram...and speech frames
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
numAntiF = 0;

plotSpecTracks2_estVar_Truth(wav, estTracks, x_errVarEKS, trueStateVTR, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileName]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

xaxe = seconds(frameInds, 1/(aParams.wLengthMS/1000*aParams.wOverlap));
xaxeS = xaxe(Sinds);
yaxeS = frameInds(Sinds);
plot(xaxeS, fs/2*yaxeS, 'g.')