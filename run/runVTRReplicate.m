% addpath(genpath('../')); % executes if directories not on path

%% Set parameters

% COMMENT OUT SAVEDIR, SOURCE, and TRACKBW WHEN RUNNING IN BATCH MODE
% output dir
savedir = '..\results\temp';

% final runs
% VTR, source 1:
% savedir = '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';
% savedir = '..\results\EKS_WS_Praat_trackBW0_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';

% VTRsynth3500, source 4:
% savedir = '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';

% VTRsynthf03500, source 5:
% savedir = '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';

% run through VTR database with given indices
vtrDbNums = 1:516; % VTRSynth
% vtrDbNums = 248:516; % EM doesn't like 247
vtrDbNums = 1;

% trackers to run 1. EKS, 2. WaveSurfer, 3. Praat
trackersToRun = 1:3;
numTrackers = length(trackersToRun);

% COMMENT OUT SOURCE WHEN RUNNING IN BATCH MODE
source = 5; % 1 - TIMIT wav files; 2 - synthesized wav files, 3 - synthesized with f0, 4 - synthesized F3 reflected, 5 - synthesized with f0 F3 reflected
cepOrder = 15;
fs = 2*3500; % if less than 8000, note that code reflects true trajectories above round(fs/2) so they are below round(fs/2)
numFormants = 3; % # formants to track
numFormants_error = 3; % # formants to do error analysis on
numAntiF = 0; % # antiformants to track
trackBW = 1; % COMMENT OUT WHEN RUNNING IN BATCH MODE
  bwFlag = 2; % applicable only if trackBW=0; 1 - Average bandwidths from WS, 2 - Fix bandwidths arbitrarily
cepData = 2; % 1 - parametric cepstrum, 2 - nonparametric cepstrum
startupframes = 0; % 20 chosen by looking at VTRsynth1 output; put back to 0 for JASA resubmission

lpcOrder = 12;   % Number of LPC Coefficients
zOrder = 0;
Dengflag = 0; % 1 - Deng linearization, 0 - Jacobian of mapping; doesn't make a big difference
useCorr = 1; % Flag for whether or not to use correlation, 1 - use, 0 - F = diag(1 1 1)
if source == 2 || source == 4 % no distinction for stochastic source
    useVAD = 0;
    VADflag = 2; % only if useVAD=1; 1 - percentile algorithm, 2 - TIMIT phones for speech/silence, 3 - TIMIT phones for phoneme classification
else
    useVAD = 1; % 0.3 s, all disp 0.6 s
    VADflag = 3; % only if useVAD=1; 1 - percentile algorithm, 2 - TIMIT phones for speech/silence, 3 - TIMIT phones for phoneme classification
end

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

rmseAll = zeros(length(vtrDbNums), numTrackers, numFormants_error+1); % database length X num algorithms (EKS, WS, Praat) X formant number (1, 2, 3, overall)
rmseAllS = zeros(length(vtrDbNums), numTrackers, numFormants_error+1); % same dims, but speech frames only
rmseAllSstartup = zeros(length(vtrDbNums), numTrackers, numFormants_error+1); % same dims, but speech frames only

if VADflag == 3 % phonetic categories
    rmseAllC = zeros(length(vtrDbNums), numTrackers, numFormants_error+1, numClasses); % similar dims, but add number of classes
    trueStateVTR2_All = cell(length(vtrDbNums), numTrackers, numFormants_error, numClasses);
    estTracks2_All = cell(length(vtrDbNums), numTrackers, numFormants_error, numClasses);
else % no phonetic categories
    trueStateVTR2_All = cell(length(vtrDbNums), numTrackers, numFormants_error);
    estTracks2_All = cell(length(vtrDbNums), numTrackers, numFormants_error);
end

% tic
for ii = vtrDbNums
    %% %%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            %wsFileName   = strcat('../data/VTR_Timit/IS07 10k resampling/Timit',num2str(vtrDbNum),'.FRM');
            wsFileName   = strcat('../data/VTR_Timit/7k resampling tracking 3 formants/Timit',num2str(vtrDbNum),'.FRM');
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
    
    % Run Snack code to find WaveSurfer tracks
    %WS_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
    params = getFormantParamsDefaultWS;
        params.numformants = numFormants;
        params.windowtype = wType;
        params.windowlength = wLengthMS/1000;
        params.framelength = params.windowlength*(1-wOverlap);
        params.lpctype = 1;
        params.lpcorder = lpcOrder;
        params.preemphasisfactor = peCoeff;
        params.ds_freq = fs;
    [trueStateWS, bwDataWS] = getFormantTrackWS(dataFileNameIn, fs_in, [], params);

%         % Read in wavesurfer .frm file
%         [trueStateWS, bwDataWS] = wavesurferFormantRead(wsFileName,numFormants);
%         %bwDataWS    = bwDataWS/2; % divide by 2 if desired to see effect
    
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
    if source >= 4 && fs < 8000
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
    % DDM 11/30/11: resample depending on WaveSurfer F3...no doesn't help
%     maxF = max(max(trueStateWS));
%     fs = 2*(100*ceil(maxF/100)+500); % round to nearest 100, then add 100 to get bandwidth
    %fs = 2*5500;
    
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
    if useVAD
        frameInds = frameInds(:, 1:numFrames)';
    else
        frameInds = formantInds(:,1);
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
    QscaleF = 5e4; %90e4, 5e4, 3e4, 1e4, 1000^2; %50e4;
    if (trackBW)
        % Get state variance from read in tracks
        %Q = diag([var(trueStateWS,0,2); var(bwDataWS,0,2)]/1);
        %Q = diag([var(diff(trueStateWS(:,frameInds>0),1,2),0,2); var(diff(bwDataWS(:,frameInds>0),1,2),0,2)]/1);

        % fix for formants and bandwidths
        QscaleBW = 1e4;
        Q = diag([QscaleF*ones(numFormants,1); QscaleBW*ones(numFormants,1)]);
        %         Q = [5e4 0 0 0 0 0;
        %             0 50e4 0 0 0 0;
        %             0 0 5e4 0 0 0;
        %             0 0 0 1e4 0 0;
        %             0 0 0 0 1e4 0;
        %             0 0 0 0 0 1e4];
    else
        % Currently is estimated from the loaded data from Wavesurfer
        %Qscale = 1; % For experiments scaling this value
        %Q = diag(var(trueState,0,2))/Qscale; %true variance PER TRACK
        %Q = eye(numFormants)*var(trueStateWS(formantInds'==1)); % to match Interspeech 2007 code in MyTrackExp.m, SAME VARIANCE FOR EACH TRACK
        %Q = diag(var(trueStateWS,0,2)); %WS variance PER TRACK
        %Q = diag(var(diff(trueStateWS(:,frameInds>0),1,2),0,2)); %WS derivative variance PER TRACK
        
        % fix for formants 
        Q = QscaleF*eye(numFormants); % set as constant across states
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
    initFormant = 500 + 1000*(0:(numFormants - 1))';
    
    % 11/1/11: initialize using WS knowledge at first speech frame
%     Sinds = find(frameInds > 0);
%     initFormant = trueStateWS(:,Sinds(1));    
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
    
    countTrack = 1;
    rmse    = zeros(numFormants_error, sum(algFlag));
    relRmse = zeros(numFormants_error, sum(algFlag));
    rmseS    = zeros(numFormants_error, sum(algFlag));
    relRmseS = zeros(numFormants_error, sum(algFlag));

    rmseSstartup    = zeros(numFormants_error, sum(algFlag));
    
    if VADflag == 3
        rmseC    = zeros(numFormants_error, numClasses);
        rmseAllC = zeros(numFormants_error, numFormants_error, numClasses);
    end
    
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
                 params.timestep = params.winLen*(1-aParams.wOverlap);
                 params.maxformantHz = fs/2;
                 params.peHz = params.maxformantHz;
                 if (source == 1 && (vtrDbNum == 153 || vtrDbNum == 193)) || ...
                     (source == 3 && (vtrDbNum == 78 || vtrDbNum == 464)) || ...
                     (source == 5 && (vtrDbNum == 124 || vtrDbNum == 131 || vtrDbNum == 138 || vtrDbNum == 217 || vtrDbNum == 228 || vtrDbNum == 416))
                     params.method = 'sl';
                 end
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
            s = 0;
            for j = 1:numFormants_error
                rmse(j,countTrack) = sqrt(mean((estTracks2(j,:,countTrack)-trueStateVTR2(j,:)).^2));

                % Compute errors where there was speech energy (i.e., omit silences)
                Sinds = find(formantInds(:,j) == 1);
                rmseS(j,countTrack) = sqrt(mean((estTracks2(j,Sinds(1:end),countTrack)-trueStateVTR2(j,Sinds(1:end))).^2)); 
                rmseSstartup(j,countTrack) = sqrt(mean((estTracks2(j,Sinds(startupframes+1:end),countTrack)-trueStateVTR2(j,Sinds(startupframes+1:end))).^2));

                rmseAll(vtrDbNum, trackerType, j) = rmse(j,countTrack);
                rmseAllS(vtrDbNum, trackerType, j) = rmseS(j,countTrack);
                rmseAllSstartup(vtrDbNum, trackerType, j) = rmseSstartup(j,countTrack);

                trueStateVTR2_All{vtrDbNum,trackerType,j} = trueStateVTR2(j,Sinds);
                estTracks2_All{vtrDbNum,trackerType,j} = estTracks2(j,Sinds);

                sDisp  = sprintf('Formant %d RMSE: %2.2f; NS RMSE: %2.2f\n', j, rmse(j,countTrack), rmseS(j,countTrack));
                s = [s sDisp];            
            end

            rmseAll(vtrDbNum, trackerType, numFormants_error+1) = mean(rmse(:,countTrack));
            rmseAllS(vtrDbNum, trackerType, numFormants_error+1) = mean(rmseS(:,countTrack));
            rmseAllSstartup(vtrDbNum, trackerType, numFormants_error+1) = mean(rmseSstartup(:,countTrack));

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
                    Cinds = find(frameInds == class); % ignore startup here, if desired; not doing it
                    rmseC(j,class) = sqrt(mean((estTracks2(j,Cinds)-trueStateVTR2(j,Cinds)).^2));
                    rmseAllC(vtrDbNum,trackerType,j,class) = rmseC(j,class,countTrack);
                    
                    trueStateVTR2_All{vtrDbNum,trackerType,j,class} = trueStateVTR2(j,Cinds);
                    estTracks2_All{vtrDbNum,trackerType,j,class} = estTracks2(j,Cinds);

                    cDisp  = sprintf('Formant %d, Class %d, RMSE: %2.2f\n', j, class, rmseC(j,class));
                    c = [c cDisp];
                end
            end
            
            if 0%verbose
                % Show individual RMSE values
                disp(c);
            end  
        end
        %^--------------^
    end
    
    if ~strcmpi(savedir,''), save(fullfile(savedir, [prefix, num2str(ii)])), end
end

%%
% break

% %% A basic plotting routine to visualize results
% plotStateTracks(trueStateVTR,estTracks,titleCell);
% 
% %% additional plotting routines Nov 2010
% figure, plot(trueStateVTR')
% hold on, plot(bwDataVTR')
% % figure, plot(trueStateWS')
% figure, plot(estTracks')
% 
% %% state tracks with uncertainties and truth
% plotStateTracksFZ_EstVar_Truth(trueStateVTR,estTracks,estVar,titleCell,numFormants,trackBW)
% disp(['Mean RMSE: ', num2str(mean(rmse))])
% rmse
% 
% %% Super-impose Ground Truth over a spectrogram
% aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
% aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
% 
% plotSpecTracks2BW(wav, trueStateVTR, aParams, numAntiF, 0);
% firstline = get(get(gca,'Title'),'String');
% title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% % set(gca, 'PlotBoxAspectRatio', [10 1 1])
% 
% %% Super-impose EKS over a spectrogram
% aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
% aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
% 
% figure
% plotSpecTracks2BW(wav, estTracks, aParams, numAntiF, trackBW);
% firstline = get(get(gca,'Title'),'String');
% title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% % set(gca, 'PlotBoxAspectRatio', [10 1 1])
% 
% %% Wavesurfer tracks on spectrogram
% figure
% plotSpecTracksWS(wav, trueStateWS, bwDataWS, aParams, trackBW)
% firstline = get(get(gca,'Title'),'String');
% title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% % set(gca, 'PlotBoxAspectRatio', [10 1 1])
%  
% %% Praat tracks on spectrogram
% figure
% praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
% params = getFormantParamsDefault;
%  params.numTracks = numFormants;
%  params.winLen = aParams.wLengthMS/1000;
%  params.timestep = params.winLen*aParams.wOverlap;
%  params.maxformantHz = fs/2;
%  params.peHz = params.maxformantHz;
% [f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(wav, aParams.fs, praat_dir, params);
% plotSpecTracksPraat(wav, f', bw', aParams, formantTimes, trackBW);
% firstline = get(get(gca,'Title'),'String');
% title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% % set(gca, 'PlotBoxAspectRatio', [10 1 1])
% 
% %% Ground truth and EKS on spectrogram
% aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
% aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
% 
% figure
% plotSpecTracks2BW_Truth(wav, trueStateVTR, estTracks, aParams, numAntiF, trackBW);
% firstline = get(get(gca,'Title'),'String');
% title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% % set(gca, 'PlotBoxAspectRatio', [10 1 1])
% 
% %% Ground truth and WaveSurfer on spectrogram
% figure
% plotSpecTracksWS_Truth(wav, trueStateWS, bwDataWS, trueStateVTR, aParams, trackBW)
% firstline = get(get(gca,'Title'),'String');
% title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% % set(gca, 'PlotBoxAspectRatio', [10 1 1])
% 
% %% Ground truth and Praat tracks on spectrogram
% praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
% params = getFormantParamsDefault;
%  params.numTracks = numFormants;
%  params.winLen = aParams.wLengthMS/1000;
%  params.timestep = params.winLen*(1-aParams.wOverlap);
%  params.maxformantHz = fs/2;
%  params.peHz = params.maxformantHz;
% [f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(wav, aParams.fs, praat_dir, params);
% plotSpecTracksPraat_Truth(wav, f', bw', trueStateVTR, aParams, formantTimes, trackBW);
% firstline = get(get(gca,'Title'),'String');
% title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% % set(gca, 'PlotBoxAspectRatio', [10 1 1])
% 
% %% EKS with uncertainty on spectrogram
% aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
% aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
% 
% figure
% plotSpecTracks2_estVar(wav, estTracks, x_errVar, aParams, numAntiF, trackBW);
% firstline = get(get(gca,'Title'),'String');
% title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% % set(gca, 'PlotBoxAspectRatio', [3 1 1])
% 
% %% Ground truth and EKS with uncertainty on spectrogram
% aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
% aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
% 
% plotSpecTracks2_estVar_Truth(wav, estTracks, x_errVar, trueStateVTR, aParams, numAntiF, trackBW);
% firstline = get(get(gca,'Title'),'String');
% title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% % set(gca, 'PlotBoxAspectRatio', [7 1 1])
% 
% %% ...and speech frames
% % Sinds(Sinds > length(frameInds)) = []; % added to enforce alignment with frameInds
% % xaxe = seconds(frameInds, 1/(aParams.wLengthMS/1000*aParams.wOverlap));
% sLength = round((numFrames+1)*wLength/2);
% xaxe = 1/fs*(wLength*(1-wOverlap):wLength*(1-wOverlap):sLength-wLength*(1-wOverlap)+1);
% xaxeS = xaxe(Sinds);
% yaxeS = frameInds(Sinds);
% plot(xaxeS, 0, 'g.')
% 
% %% Ground truth and EKS on spectrogram...and speech frames
% aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
% aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even
% 
% figure
% plotSpecTracks2BW_Truth(wav, trueStateVTR, estTracks, aParams, numAntiF, trackBW);
% % plotSpecTracks2BW(wav, estTracks, aParams, numAntiF, trackBW);
% firstline = get(get(gca,'Title'),'String');
% title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% 
% Sinds(Sinds > length(frameInds)) = []; % added to enforce alignment with frameInds
% sLength = round((numFrames+1)*wLength/2);
% xaxe = 1/fs*(wLength*(1-wOverlap):wLength*(1-wOverlap):sLength-wLength*(1-wOverlap)+1);
% xaxeS = xaxe(Sinds);
% yaxeS = frameInds(Sinds);
% plot(xaxeS, 0, 'g.')
% % ylim([0 5000])