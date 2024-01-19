% addpath(genpath('../')); % executes if directories not on path

%% Set parameters

% output dir
savedir = '..\results\temp';
% savedir = '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs10000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWws';

% run through VTR database with given indices
% vtrDbNums = [1:152 154:192 194:516]; % when doing praat on VTR, Praat bombs on 153 and 193
% vtrDbNums = 1:516; % VTRSynth
% vtrDbNums = [1:77 79:463 465:516]; % for praat on VTRsynthf0
% vtrDbNums = 248:516; % EM doesn't like 247
% vtrDbNums = [1:123 125:130 132:216 218:415 417:516]; % for praat on VTRsynthf03500
vtrDbNums = 1;
    
% trackers to run 1. EKS, 2. WaveSurfer, 3. Praat
trackersToRun = 1:3;
numTrackers = length(trackersToRun);

source = 1; % 1 - TIMIT wav files; 2 - synthesized wav files, 3 - synthesized with f0, 4 - synthesized F3 reflected, 5 - synthesized with f0 F3 reflected
cepOrder = 15;
fs = 2*3500; % if less than 8000, note that code reflects true trajectories above round(fs/2) so they are below round(fs/2)
numFormants = 3; % # formants to track
numAntiF = 0; % # antiformants to track
trackBW = 1;
  bwFlag = 2; % applicable only if trackBW=0; 1 - Average bandwidths from WS, 2 - Fix bandwidths arbitrarily
cepData = 1; % 1 - parametric cepstrum, 2 - nonparametric cepstrum
startupframes = 0; % 20 chosen by looking at VTRsynth1 output

lpcOrder   = 12;   % Number of LPC Coefficients
zOrder = 0;
Dengflag = 0; % 1 - Deng linearization, 0 - Jacobian of mapping; doesn't make a big difference
useCorr = 1; % Flag for whether or not to use correlation, 1 - use, 0 - F = diag(1 1 1)
useVAD = 1; % 0.3 s, all disp 0.6 s
  VADflag = 3; % only if useVAD=1; 1 - percentile algorithm, 2 - TIMIT phones for speech/silence, 3 - TIMIT phones for phoneme classification
  numClasses = 6; % number of phoneme classes, Deng06 has 6
  
% Select which tracker to run (pick ONE)
algFlag = [0 1 0 0]; % Select 1 to run, 0 not to
EKF = 1; EKS = 2; EKS_EM = 3; PF = 4;

if ~exist(savedir, 'dir'), mkdir(savedir), end

% prefix to filename in each directory
switch source
    case 1
        prefix = 'VTR';
        peCoeff    = 0.7;  % Pre-emphasis factor: 0 for stochastic sources, 0.7 if voicing exists
    case 2
        prefix = 'VTRsynth';
        peCoeff    = 0;  % Pre-emphasis factor: 0 for stochastic sources, 0.7 if voicing exists
    case 3
        prefix = 'VTRsynthf0';
        peCoeff    = 0.7;  % Pre-emphasis factor: 0 for stochastic sources, 0.7 if voicing exists
    case 4
        prefix = 'VTRsynth3500';
        peCoeff    = 0;  % Pre-emphasis factor: 0 for stochastic sources, 0.7 if voicing exists
    case 5
        prefix = 'VTRsynthf03500';
        peCoeff    = 0.7;  % Pre-emphasis factor: 0 for stochastic sources, 0.7 if voicing exists
end

verbose = 1; % display RMSE

% Analysis parameters
wType = 'hamming'; % Window type
wLengthMS  = 20;   % Length of window (in milliseconds)
wOverlap   = 0.5;  % Factor of overlap of window

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

rmseAll = zeros(length(vtrDbNums), numTrackers, numFormants+1); % database length X num algorithms (EKS, WS, Praat) X formant number (1, 2, 3, overall)
rmseAllS = zeros(length(vtrDbNums), numTrackers, numFormants+1); % same dims, but speech frames only
rmseAllC = zeros(length(vtrDbNums), numTrackers, numFormants+1, numClasses); % similar dims, but add number of classes

trueStateVTR2_All = cell(length(vtrDbNums), numTrackers, numFormants, numClasses);
estTracks2_All = cell(length(vtrDbNums), numTrackers, numFormants, numClasses);

%%
% tic
for ii = vtrDbNums
    %%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vtrDbNum = ii;
    %disp('Database number')
    disp(vtrDbNum)

    switch source
        case 1
            %%%%%%%%%%%%% With TIMIT sentences %%%%%%%%%%%%%%%%%%
            % Waveform path
            dataFileName = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.wav');
            dataFileNameIn = dataFileName;

            % Read in audio waveform
            [wav_orig, fs_in] = wavread(dataFileName);

            % Wavesurfer file path
            %wsFileName   = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.FRM');
            %wsFileName   = strcat('../data/VTR_Timit/tcl params 4 Hamming .02 .01 1 12 .7 10000 500/Timit',num2str(vtrDbNum),'.FRM');
            wsFileName   = strcat('../data/VTR_Timit/IS07 10k resampling/Timit',num2str(vtrDbNum),'.FRM');
            %wsFileName   = strcat('../data/VTR_Timit/7k resampling tracking 3 formants/Timit',num2str(vtrDbNum),'.FRM');
        case 2
            %%%%%%%%% With synthesized waveforms using VTR trajectories %%%%
            % Waveform path
            dataFileName = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.wav');
            dataFileNameIn = strcat('../data/VTRsynth/VTRsynth',num2str(vtrDbNum),'.wav');

            % Read in audio waveform
            [wav_orig, fs_in] = wavread(dataFileNameIn);

            % Wavesurfer file path
            wsFileName   = strcat('../data/VTRsynth/tcl params 4 Hamming .02 .01 1 12 0 10000 500/VTRsynth',num2str(vtrDbNum),'.FRM');
        case 3
            %%%%%%%%% With synthesized waveforms using VTR trajectories and f0 %%%%
            % Waveform path
            dataFileName = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.wav');
            dataFileNameIn = strcat('../data/VTRsynthf0/VTRsynth',num2str(vtrDbNum),'.wav');

            % Read in audio waveform
            [wav_orig, fs_in] = wavread(dataFileNameIn);

            % Wavesurfer file path
            %wsFileName   = strcat('../data/VTRsynthf0/tcl params 4 Hamming .02 .01 1 12 0.7 16000 500/VTRsynth',num2str(vtrDbNum),'.FRM');
            wsFileName   = strcat('../data/VTRsynthf0/tcl params 3 Hamming .02 .01 1 12 0.7 7000 500/VTRsynth',num2str(vtrDbNum),'.FRM');
        case 4
            %%%%%%%%% With synthesized waveforms using VTR trajectories, F3 reflected under 3500 Hz %%%%
            % Waveform path
            dataFileName = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.wav');
            dataFileNameIn = strcat('../data/VTRsynth - F3 under 3500/VTRsynth',num2str(vtrDbNum),'.wav');

            % Read in audio waveform
            [wav_orig, fs_in] = wavread(dataFileNameIn);

            % Wavesurfer file path
            wsFileName   = strcat('../data/VTRsynth - F3 under 3500/tcl params 3 Hamming .02 .01 1 12 0 7000 500/VTRsynth',num2str(vtrDbNum),'.FRM');
        case 5
            %%%%%%%%% With synthesized waveforms using VTR trajectories and f0, F3 reflected under 3500 Hz %%%%
            % Waveform path
            dataFileName = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.wav');
            dataFileNameIn = strcat('../data/VTRsynthf0 - F3 under 3500/VTRsynth',num2str(vtrDbNum),'.wav');

            % Read in audio waveform
            [wav_orig, fs_in] = wavread(dataFileNameIn);

            % Wavesurfer file path
            wsFileName   = strcat('../data/VTRsynthf0 - F3 under 3500/tcl params 3 Hamming .02 .01 1 12 0.7 7000 500/VTRsynth',num2str(vtrDbNum),'.FRM');            
        otherwise
            error('Incorrect source.')
    end
    
    % Read in wavesurfer waveform
    [trueStateWS, bwDataWS] = wavesurferFormantRead(wsFileName,numFormants);
    %bwDataWS    = bwDataWS/2; % divide by 2 if desired to see effect
    clear data
    
    % Read in 'Truth' values from VTR database
    load ../data/VTR_Timit/allVTRTracks;
    curData = DATA{vtrDbNum};
    clear DATA

    % Rescale appropriately
    trueStateVTR = 1000*curData.vtrData(:,1:numFormants)';    % Scaling from units
    
    if source > 1
        % truth for VTRsynth and VTRsynthf0 (what I synthesized with)
        bwDataVTR    = 500*curData.vtrData(:,5:(4+numFormants))';
    else
        % truth for VTR database, given in kHz, so need to multiply by 1000
        bwDataVTR    = 1000*curData.vtrData(:,5:(4+numFormants))';
    end
    
    % 6/20/11 modification to reflect F3 under 3500 Hz, only for
    % synthesized waveforms
    if source > 1 && fs < 8000
        cutoff = round(fs/2);
        F3 = trueStateVTR(3, 1:end);
        [temp, clipindex] = find(F3 > cutoff);
        F3clip = F3(clipindex);
        F3noclip = cutoff - (F3clip - cutoff);
        F3(clipindex) = F3noclip;
        trueStateVTR(3, 1:end) = F3; % modify formant matrix
        clear temp
    end

    %%%%%%%%%%%%%%%%%%%%%% Generate observation sequence %%%%%%%%%%%%%%%%%%%%%%
    % Resample input data if sampling rate is not equal to that of input
    if(fs ~= fs_in)
        wav = resample(wav_orig,fs,fs_in,2048);
        %display(['Input Fs = ' int2str(fs_in) ' Hz; resampling to ' num2str(fs) ' Hz'])
    else
        wav = wav_orig;
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
    switch cepData
        case 1
            % y = genLPCC(wav, win, wOverlap, peCoeff, lpcOrder, cepOrder);
            y = genLPCCz(wav, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);
        case 2
            y = genCeps(wav, win, wOverlap, peCoeff, cepOrder); % nonparametric cepstrum
        otherwise
            error('Error in switch cepData.')
    end

    % Truncate at most one frame from output of LPCC call
    numFrames    = min([size(y,2), size(trueStateWS, 2), size(trueStateVTR, 2)]);
    numObs       = numFrames;
    
    y            = y(:, 1:numFrames);
    trueStateWS  = trueStateWS(1:numFormants,1:numFrames);
    bwDataWS     = bwDataWS(1:numFormants,1:numFrames);
    trueStateVTR = trueStateVTR(1:numFormants,1:numFrames);
    bwDataVTR    = bwDataVTR(1:numFormants,1:numFrames);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Speech Activity Detection/Classification %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(useVAD)
        %display('Using Voice Activity Detection');
        % If all formants are coasted jointly, then do not do multiband detection
        multiBand = ~coastJoint;
        
        switch VADflag
            case 1
                % simple VAD
                frameInds = simpleVAD(wav, aParams, quantThresh, multiBand, [], plotVad);
                %frameInds = zeros(1, 247);
                %frameInds(13:154) = 1;                
            case 2
                % TIMIT labels to help define speech frames
                load([dataFileName(1:end-4), '.mat']) % 0.44 s
                frameInds = timitVAD(wav, data.phones, aParams, plotVad); % 0.3 s
            case 3
                load([dataFileName(1:end-4), '.mat']) % 0.44 s
                frameInds = TIMITphnclassify(wav, data.phones, aParams, plotVad); % 0.3 s
            otherwise
                error('Incorrect VADflag selected.')
        end
        
        if(multiBand)
            if(numFormants > 4)
                error('Multi-band VAD not supported for more than 4 formants');
            else
                % Formant energy bands (taken from Deng et al, 2006)
                bandEdges = [200 900; 600 2800; 1400 3800; 1700 5000];
                bandEdges = bandEdges(1:numFormants,:);
            end

            formantInds = logical(frameInds);
        else
            if trackBW
                formantInds = repmat(logical(frameInds),2*numFormants,1)';
            else
                formantInds = repmat(logical(frameInds),numFormants,1)';
            end
        end
    else
        if trackBW
            formantInds = ones(numFrames,2*numFormants);
        else
            formantInds = ones(numFrames,numFormants);
        end
    end

    % get rid of indices that are out of range
    formantInds = formantInds(1:numFrames, :);
    frameInds = frameInds(:, 1:numFrames)';
    
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
        %Q = diag([var(trueStateWS,0,2); var(bwDataWS,0,2)]/1);
        %Q = diag([var(diff(trueStateWS(:,frameInds>0),1,2),0,2); var(diff(bwDataWS(:,frameInds>0),1,2),0,2)]/1);
        
        % fix for formants and bandwidths
        QscaleF = 10e4;
        QscaleBW = 1e4;
        Q = diag([QscaleF*ones(numFormants,1); QscaleBW*ones(numFormants,1)]);
    else
        % Currently is estimated from the loaded data from Wavesurfer
        %Qscale = 1; % For experiments scaling this value
        %Q = diag(var(trueState,0,2))/Qscale; %true variance PER TRACK
        %Q = eye(numFormants)*var(trueStateWS(formantInds'==1)); % to match Interspeech 2007 code in MyTrackExp.m, SAME VARIANCE FOR EACH TRACK
        %Q = diag(var(trueStateWS,0,2)); %WS variance PER TRACK
        Q = diag(var(diff(trueStateWS(:,frameInds>0),1,2),0,2)); %WS derivative variance PER TRACK
        
        %Qscale = QscaleV(rr);
        
        % Q for each state derived from WS tracks
        %Qscale = 5e4;
        %Q = Qscale*eye(numFormants); % set as constant across states
    end

    %%% Set measurement noise covariance matrix R %%%%

    % Observation noise ``variance'', should decrease as the cepstral order increases
    sigExp = 1;
    R      = diag(1./ones(cepOrder,1)./(((1:cepOrder).^sigExp)'));
    %R = 1/10*diag(1./ones(cepOrder,1));

    %%% Set process evolution matrix F %%%
    % Estimate cross-correlation among Formant Trajectories
    if(useCorr)
        if(diagCorr)
            trueStateWStmp = trueStateWS(:,formantInds(:,1)==1);
            for dim = 1:numFormants
                ar(dim) = fitVAR(trueStateWStmp(dim,:)',1);
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
        switch bwFlag 
            case 1 % fix bandwidths from averaged WaveSurfer tracks
                bwDataWStmp = bwDataWS(:,formantInds(:,1)==1);
                bwStates = repmat(mean(bwDataWStmp,2),1,size(bwDataWS,2));
            case 2 % fix bandwidths from arbitrary settings
                fixBW = 80 + 40*(0:(numFormants - 1))';
                bwStates = repmat(fixBW,1,numFrames);
            otherwise
                error('Incorrect bwFlag setting.')
        end
    end

    %%% Set initial state of formant trackers
    % Using values suggested by Li Deng and otherwise in literature
    %initFormant = 500 + 1000*(0:(numFormants - 1))';
    
    % 11/1/11: initialize using WS knowledge at first speech frame
    Sinds = find(frameInds > 0);
    initFormant = trueStateWS(:,Sinds(1));    
    initBW      = 80 + 40*(0:(numFormants - 1))';
    if(trackBW)
        x0 = [initFormant', initBW']';
    else
        x0 = initFormant;
    end

    %%% TRACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    countTrack = 1; % Counter for storing results
    % Initialize root-mean-square error matrices:
    rmse    = zeros(numFormants, sum(algFlag));
    relRmse = zeros(numFormants, sum(algFlag));
    clear estTracks estVar

    % Run Extended Kalman Filter
    if algFlag(EKF)
        smooth = 0;  % No smoothing
        %[x_est x_errVar] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth);
        [x_est x_errVar] = formantTrackEKSZ(y, F, Q, R, x0, formantInds', fs, bwStates, numFormants, smooth, Dengflag);

        %Track estimate into data cube for plot routines
        estTracks(:,:,countTrack) = x_est;
        estVar(:,:,:,countTrack) = x_errVar;
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
        %[x_est x_errVar] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth);
        %[x_est x_errVar] = formantTrackEKSZ(y, F, Q, R, x0, formantInds', fs, bwStates, numFormants, smooth);
        [x_est x_errVar] = formantTrackEKSZ(y, F, Q, R, x0, formantInds', fs, bwStates, numFormants, smooth, Dengflag);

        % Track estimate into data cube for plot routines
        estTracks(:,:,countTrack) = x_est;
        estVar(:,:,:,countTrack) = x_errVar;
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

%         % Compute and display RMSE and relative RMSE
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
%         countTrack = countTrack + 1;    % Increment counter
    end


    if algFlag(PF)
        numParticles = 1000;
        [x_est x_errVar] = formantTrackPF(y, F, Q, R, x0, formantInds, fs, bwStates, numParticles);

        % Track estimate into data cube for plot routines
        estTracks(:,:,countTrack) = x_est;
        estVar(:,:,:,countTrack) = x_errVar;
        titleCell(1,countTrack+1) = {'PF'};
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
%         disp('Particle Filter Results');
%         disp(s); % Show individual RMSE values
%         disp(sprintf('Average RMSE: %2.2f; NS RMSE: %2.2f', mean(rmse(:,countTrack)),mean(rmseS(:,countTrack))));
%         countTrack = countTrack + 1;
    end


    %Initial Plotting Variables
    titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
    titleCell(2,1)  = {'r'};            % Color for true state plot
    
    %% find RMSE, run this cell after loading a given .mat file
    %% (e.g., VTR200.mat) with verbose = 1, but RMSE determined within
    %% each formant and then average RMSE found
    countTrack = 1;
    numFormants_error = 3;
    rmse    = zeros(numFormants_error, sum(algFlag));
    relRmse = zeros(numFormants_error, sum(algFlag));
    rmseS    = zeros(numFormants_error, sum(algFlag));
    relRmseS = zeros(numFormants_error, sum(algFlag));

    rmseC    = zeros(numFormants_error, numClasses);
    rmseAllC = zeros(numFormants_error, numFormants_error, numClasses);
    
    % choose tracker: 1. EKS, 2. WaveSurfer, 3. Praat
    for jj = trackersToRun
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
                frameInds = frameInds(1:len, :);
                algo = 'Praat Results';
            otherwise
                error('Unknown tracker type.')
        end

        %Compute and Display MSE and RMSE
        %trueStateVTR2 = trueStateVTR;
        s = 0;
        for j = 1:numFormants_error
            rmse(j,countTrack) = sqrt(mean((estTracks2(j,:,countTrack)-trueStateVTR2(j,:)).^2));

            % Compute errors where there was speech energy (i.e., omit silences)
            Sinds = find(formantInds(:,j) == 1);
            rmseS(j,countTrack) = sqrt(mean((estTracks2(j,Sinds(startupframes+1:end),countTrack)-trueStateVTR2(j,Sinds(startupframes+1:end))).^2)); % ignore startup

            rmseAll(vtrDbNum, trackerType, j) = rmse(j,countTrack);
            rmseAllS(vtrDbNum, trackerType, j) = rmseS(j,countTrack);
            
            sDisp  = sprintf('Formant %d RMSE: %2.2f; NS RMSE: %2.2f\n', j, rmse(j,countTrack), rmseS(j,countTrack));
            s = [s sDisp];            
        end
        rmseAll(vtrDbNum, trackerType, numFormants_error+1) = mean(rmse(:,countTrack));
        rmseAllS(vtrDbNum, trackerType, numFormants_error+1) = mean(rmseS(:,countTrack));
        
        sDisp  = sprintf('Average   RMSE: %2.2f; NS RMSE: %2.2f\n', mean(rmse(:,countTrack)),mean(rmseS(:,countTrack)));
        s = [s sDisp];
            
        if verbose
            % Show individual RMSE values
            disp(algo);
            disp(s)
        end
 
        %v------Now find RMSE per phoneme class (6 classes plus silence label possible) --------v        
        if VADflag == 3
            c = 0;
            for j = 1:numFormants_error
                % Compute errors for particular phoneme classes
                % class = silence (0), vowel (1), semivowel (2), nasal (3), fricative (4),
                % affricate (5), stop (6)
                for class = 1:numClasses
                    Cinds = find(frameInds == class);
                    rmseC(j,class) = sqrt(mean((estTracks2(j,Cinds(startupframes+1:end))-trueStateVTR2(j,Cinds(startupframes+1:end))).^2)); % ignore startup
                    rmseAllC(vtrDbNum,trackerType,j,class) = rmseC(j,class,countTrack);
                    
                    trueStateVTR2_All{vtrDbNum,trackerType,j,class} = trueStateVTR2(j,Cinds);
                    estTracks2_All{vtrDbNum,trackerType,j,class} = estTracks2(j,Cinds);

                    cDisp  = sprintf('Formant %d, Class %d, RMSE: %2.2f\n', j, class, rmseC(j,class));
                    c = [c cDisp];
                end
            end
            
            if verbose
                % Show individual RMSE values
                disp(c);
            end
        end
        %^--------------^
        
    end
    
    %if ~strcmpi(savedir,''), save(fullfile(savedir, [prefix, num2str(ii)])), end
end

%%
break
%% Calculate error, used for JASA submission outputs
% load('..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWws\VTR516.mat')
% indices = [1:152 154:192 194:516]; % for VTR (praat bombs)
% indices = 1:516; % for VTRSynth
% indices = [1:77 79:463 465:516]; % for praat
indices = [1:123 125:130 132:216 218:415 417:516]; % for praat on VTRsynthf03500

% rmseAllS(vtrDbNum, trackerType, numFormants_error+1)
% 1. EKS, 2. WaveSurfer, 3. Praat
% indices = 1:516;
rmse_table = zeros(numTrackers, numFormants_error+1);
for ii = trackersToRun % # trackers
    for jj = 1:4 % # rmse values
        tmp = rmseAllS(indices, ii, jj);
        rmse_table(ii, jj) = mean(tmp(~isnan(tmp))); % take into account any NaNs
    end
end

figure, plot(indices, rmseAllS(indices, 1, 4), 'o')

rmse_table'
round(rmse_table')

%% For JASA re-submission: Calculate error within phonetic classes
% load('..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWws\VTR516.mat')

% vtrDbNums = [1:152 154:192 194:516]; % when doing praat, Praat bombs on 153 and 193

trueStateVTR2_All_Db = cell(numTrackers, numFormants_error, numClasses);
estTracks2_All_Db = cell(numTrackers, numFormants_error, numClasses);
trueStateVTR2_All_Db_class = cell(numTrackers, numFormants_error);
estTracks2_All_Db_class = cell(numTrackers, numFormants_error);
trueStateVTR2_All_Db_class_f = cell(1, numTrackers);
estTracks2_All_Db_class_f = cell(1, numTrackers);
for trackerType = trackersToRun
    for j = 1:numFormants_error
        for class = 1:numClasses
            for vtrDbNum = vtrDbNum % *** set to vtrDbNum if only interested in current utterance RMSE
                % reshape so that we collapse tracks across all utterances, preserving
                % categories of formant number and phoneme class (and
                % tracker type)
                trueStateVTR2_All_Db{trackerType,j,class} = [trueStateVTR2_All_Db{trackerType,j,class} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                estTracks2_All_Db{trackerType,j,class}    = [estTracks2_All_Db{trackerType,j,class}    estTracks2_All{vtrDbNum,trackerType,j,class}];
                
                % reshape so that we collapse tracks across all utterances and all phonetic classes, preserving
                % category of formant number (and tracker type)                
                trueStateVTR2_All_Db_class{trackerType,j} = [trueStateVTR2_All_Db_class{trackerType,j} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                estTracks2_All_Db_class{trackerType,j}    = [estTracks2_All_Db_class{trackerType,j}    estTracks2_All{vtrDbNum,trackerType,j,class}];                
                
                % reshape so that we collapse tracks across all utterances, all phonetic classes, and all formants, preserving
                % tracker type
                trueStateVTR2_All_Db_class_f{trackerType} = [trueStateVTR2_All_Db_class_f{trackerType} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                estTracks2_All_Db_class_f{trackerType}    = [estTracks2_All_Db_class_f{trackerType}    estTracks2_All{vtrDbNum,trackerType,j,class}];                
            end
        end
    end
end

%% Run this cell after loading .mat file (e.g., VTR516.mat or
%% VTRsynthf0516.mat) and running previous cell
% Number of frames for each phonetic class
for ii = 1:6, temp = trueStateVTR2_All_Db{1,1,ii};disp(length(temp)), end

% Number of frames for each formant (should be the same)
for ii = 1:3, temp = trueStateVTR2_All_Db_class{1,ii};disp(length(temp)), end

% Number of frames for each tracker
temp = trueStateVTR2_All_Db_class_f;disp(length(temp))

% now, find error for each formant trajectory and put in table a la Deng (2006)
% Row 1 (Class 1): F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
% Row 2 (Class 2): F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
% ...
rmse_table = zeros(numClasses, numTrackers*numFormants_error);
for class = 1:numClasses
    for trackerType = trackersToRun
        for j = 1:numFormants_error
            differ = estTracks2_All_Db{trackerType,j,class} - trueStateVTR2_All_Db{trackerType,j,class};
            differ = differ(~isnan(differ)); % get rid of NaNs
            % calculate RMSE
            rmse_table(class, j+(trackerType-1)*numTrackers) = sqrt(mean((differ).^2));
            
            % calculate mean absolute error
            %rmse_table(class, j+(trackerType-1)*numTrackers) = mean(abs(differ));
        end
    end
end
round(rmse_table)

% now, find error for each formant trajectory and put in table
% F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
rmse_table2 = zeros(1, numTrackers*numFormants_error);
for trackerType = trackersToRun
    for j = 1:numFormants_error
        differ = estTracks2_All_Db_class{trackerType,j} - trueStateVTR2_All_Db_class{trackerType,j};
        differ = differ(~isnan(differ)); % get rid of NaNs
        % calculate RMSE
        rmse_table2(j+(trackerType-1)*numTrackers) = sqrt(mean((differ).^2));

        % calculate mean absolute error
        %rmse_table2(j+(trackerType-1)*numTrackers) = mean(abs(differ));
    end
end
round(rmse_table2)

% now, find overall error for each tracker and put in table
% F1, F2, F3 (for all trackers)
rmse_table3 = zeros(1, numTrackers);
for trackerType = trackersToRun
    differ = estTracks2_All_Db_class_f{trackerType} - trueStateVTR2_All_Db_class_f{trackerType};
    differ = differ(~isnan(differ)); % get rid of NaNs
    % calculate RMSE
    rmse_table3(trackerType) = sqrt(mean((differ).^2));

    % calculate mean absolute error
    %rmse_table3(trackerType) = mean(abs(differ));
end
round(rmse_table3)

% save(fullfile(savedir, [prefix, '_error']))

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

plotSpecTracks2BW(wav, trueStateVTR, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileNameIn]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Super-impose EKS over a spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2BW(wav, estTracks, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileNameIn]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Wavesurfer tracks on spectrogram
% [f,bw] = wavesurferFormantRead([wsFileName(1:end-3) 'frm'],numFormants);
% plotSpecTracksWS(wav, f, bw, aParams, trackBW)
plotSpecTracksWS(wav, trueStateWS, bwDataWS, aParams, trackBW)
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileNameIn]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])
 
%% Praat tracks on spectrogram
praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
params = getFormantParamsDefault;
 params.numTracks = numFormants;
 params.winLen = aParams.wLengthMS/1000;
 params.timestep = params.winLen*aParams.wOverlap;
 params.maxformantHz = fs/2;
 params.peHz = params.maxformantHz;
[f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(wav, aParams.fs, praat_dir, params);
plotSpecTracksPraat(wav, f', bw', aParams, formantTimes, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileNameIn]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and EKS on spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2BW_Truth(wav, trueStateVTR, estTracks, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileNameIn]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and WaveSurfer on spectrogram
[f,bw] = wavesurferFormantRead([wsFileName(1:end-3) 'frm'],numFormants);
plotSpecTracksWS_Truth(wav, f, bw, trueStateVTR, aParams, trackBW)
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileNameIn]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and Praat tracks on spectrogram
praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
params = getFormantParamsDefault;
 params.numTracks = numFormants;
 params.winLen = aParams.wLengthMS/1000;
 params.timestep = params.winLen*aParams.wOverlap;
 params.maxformantHz = fs/2;
 params.peHz = params.maxformantHz;
[f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(wav, aParams.fs, praat_dir, params);
plotSpecTracksPraat_Truth(wav, f', bw', trueStateVTR, aParams, formantTimes, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileNameIn]; firstline})
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% EKS with uncertainty on spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2_estVar(wav, estTracks, x_errVar, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileNameIn]; firstline})
% set(gca, 'PlotBoxAspectRatio', [3 1 1])

%% Ground truth and EKS with uncertainty on spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2_estVar_Truth(wav, estTracks, x_errVar, trueStateVTR, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileNameIn]; firstline})
% set(gca, 'PlotBoxAspectRatio', [7 1 1])

%% ...and speech frames
xaxe = seconds(frameInds, 1/(aParams.wLengthMS/1000*aParams.wOverlap));
xaxeS = xaxe(Sinds);
yaxeS = frameInds(Sinds);
plot(xaxeS, 0, 'g.')

%% Ground truth and EKS on spectrogram...and speech frames
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2BW_Truth(wav, trueStateVTR, estTracks, aParams, numAntiF, trackBW);
% plotSpecTracks2BW(wav, estTracks, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['File: ', dataFileNameIn]; firstline})

xaxe = seconds(frameInds, 1/(aParams.wLengthMS/1000*aParams.wOverlap));
xaxeS = xaxe(Sinds);
yaxeS = frameInds(Sinds);
plot(xaxeS, 0, 'g.')
ylim([0 5000])