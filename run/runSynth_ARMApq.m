function varargout = runSynth_ARMApq(F, Fbw, Z, Zbw, dur, pNoiseVar, snr_dB, cepOrder, fs, trackBW, plot_flag, algFlag, x0, aParams)

% Track poles and zeros (no bandwidths) on synthetic data
% Based on runSynthZ(), which was model-based. This function is
% waveform-based.
% 
% Author: Daryush Mehta
% Created:  04/30/2010
% Modified: 05/01/2010 variable output number
%           05/09/2010 track bandwidths
%           06/12/2010 P matrix output
%           06/16/2010 (analysis parameters on input, x output)
% 
% INPUT:
%    F:         center frequencies of the resonances (col vector), in Hz
%    Fbw:       corresponding bandwidths of the resonances (col vector), in Hz
%    Z:         center frequencies of the anti-resonances (col vector), in Hz
%    Zbw:       corresponding bandwidths of the anti-resonances (col vector), in Hz
%    dur:       duration of signal, in s
%    pNoiseVar: process noise variance
%    snr_dB:    observation noise, in dB
%    cepOrder:  Number of cepstal coefficients to compute
%    fs:        sampling rate of waveform, in Hz
%    trackBW:   track bandwidths if 1
%    plot_flag: plot figures if 1
%    algFlag:   select 1 to run, 0 not to for [EKF EKS]
%    x0:        initial state of formant trackers [F; Fbw; Z; Zbw], in Hz
%    aParams:   structure with following fields:
%                     wType       % window type
%                     wLengthMS   % Length of window (in milliseconds)
%                     wOverlap    % Factor of overlap of window
%                     lpcOrder    % Number of AR coefficients
%                     zOrder      % Number of MA coefficients
%                     peCoeff     % Pre-emphasis factor
% 
% OUTPUT:
%    Depends on algFlag. For each algorithm, three outputs then waveform--
%       rmse_mean:  average RMSE across all tracks
%       x_est:      estimated tracks
%       x_errVar:   covariance matrix of estimated tracks
%       So that if two algorithms run, the following are output:
%       [rmse_meanEKF, x_estEKF, x_errVarEKF, rmse_meanEKS, x_estEKS, x_errVarEKS, x]
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% Synthetic Examples:
%   rmse_mean = runSynth_ARMApq([500; 1000],[100; 100],[],[],.5,200,25,15,16000,1,[1 0],[600; 1600]);
%
% The usual examples from PRAAT and Wavesurfer are not included, because
% it is not clear how to generate paths of zeros along with the paths of
% formants that are available.

addpath(genpath('../')); % Paths
rand('state',sum(100*clock)); randn('state',sum(100*clock)); % Seeds

%% Create an ARMA model by filtering a white noise sequence
N = round(dur*fs);
[num, denom] = fb2tf(F, Fbw, Z, Zbw, fs);
x = filter(num, denom, randn(N,1));

if plot_flag
    % Compute transfer function coefficients
    figure

    subplot(311)
    [spec_tf, freq] = freqz(1, denom, 512, fs);
    spec_tf = 20*log10(abs(spec_tf));
    plot(freq, spec_tf)
    title('AR-only spectrum')

    subplot(312)
    [spec_tf, freq] = freqz(num, 1, 512, fs);
    spec_tf = 20*log10(abs(spec_tf));
    plot(freq, spec_tf)
    title('MA-only spectrum')

    subplot(313)
    [spec_tf, freq] = freqz(num, denom, 512, fs);
    spec_tf = 20*log10(abs(spec_tf));
    plot(freq, spec_tf)
    title('Transfer function')
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')

    %%
    figure, subplot(211)
    plot(x)
    title('Time-domain waveform');
    xlabel('Samples'); ylabel('Amplitude');

    subplot(212), hold on
    [spec_true, freq_true] = pwelch(x-mean(x), [], [], [], fs);
    spec_true = 10*log10(spec_true);
    plot(freq_true, spec_true)
    plot(freq, spec_tf, 'r')
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    legend('Periodogram', 'Transfer function')

    %%
    if ~isempty(F)
        % Estimate AR params using arcov (this leads to a biased estimate)
        [arCoeffs e] = arcov(x,length(F)*2);
        disp(' ')
        disp('True Coefficients');
        disp(['AR Coeffs:' num2str(denom)])
        disp(['MA Coeffs:' num2str(num)])
        disp('Covariance Method: Estimated AR Coefficients')
        disp(['AR Coeffs:' num2str(arCoeffs)])
        [spec, freq] = freqz(sqrt(e), arCoeffs, 512, fs);
        figure, subplot(311), hold on
        plot(freq, 20*log10(abs(spec)), 'b')
        [spec, freq] = freqz(1, denom, 512, fs);
        plot(freq, 20*log10(abs(spec)), 'r:')
        title('AR estimate spectrum (ARCOV)')
        legend('Estimate', 'True AR')
    end
    
    %% Estimate ARMA parameters using armax function from Sys. ID. toolbox
    data = iddata(x,[],1); % Package input
    m = armax(data,[length(F)*2 length(Z)*2]); % Call estimator with desired model orders
    disp('Sys ID toolbox ARMA estimates');
    disp(['AR Coeffs:' num2str(m.a)]); % Estimated AR Coefficients
    disp(['MA Coeffs:' num2str(m.c)]); % Estimated MA Coefficients

    figure

    subplot(311), hold on
    [spec, freq] = freqz(1, m.a, 512, fs);
    plot(freq, 20*log10(abs(spec)), 'b')
    [spec, freq] = freqz(1, denom, 512, fs);
    plot(freq, 20*log10(abs(spec)), 'r:')
    title('AR estimate spectrum (ARMA)')
    legend('Estimate', 'True')

    subplot(312), hold on
    [spec, freq] = freqz(m.c, 1, 512, fs);
    plot(freq, 20*log10(abs(spec)), 'b')
    [spec, freq] = freqz(num, 1, 512, fs);
    plot(freq, 20*log10(abs(spec)), 'r:')
    title('MA estimate spectrum (ARMA)')
    legend('Estimate', 'True')

    subplot(313), hold on
    [spec, freq] = freqz(m.c, m.a, 512, fs);
    plot(freq, 20*log10(abs(spec)), 'b')
    [spec, freq] = freqz(num, denom, 512, fs);
    plot(freq, 20*log10(abs(spec)), 'r:')
    title('ARMA estimate spectrum (ARMA)')
    legend('Estimate', 'True')
end

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

% uncomment next line to input truth
% y = genLPCCz(x, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder, num, denom);
c = genLPCCz(x, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);

% observation noise added given SNR in input
[y, oNoiseVar] = addONoise(c, snr_dB);

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
nP = length(F);
nZ = length(Z);
numObs = size(y, 2);

if trackBW
    trueState = [repmat(F, 1, numObs); repmat(Fbw, 1, numObs); ...
        repmat(Z, 1, numObs); repmat(Zbw, 1, numObs)];
    bwStates = []; % If we are tracking bandwidths do not provide them
else
    trueState = repmat([F; Z], 1, numObs);
    bwStates = repmat([Fbw; Zbw], 1, numObs);
end
stateLen = size(trueState,1);

% Process Matrix F
Fmatrix = eye(stateLen);

% process noise is input parameter
Q = diag(pNoiseVar*ones(stateLen,1));

R = oNoiseVar*eye(cepOrder); % Measurement noise covariance matrix R perfectly matches added noise

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
    varargout(countOut) = {rmse_mean}; countOut = countOut + 1;
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
    varargout(countOut) = {rmse_mean}; countOut = countOut + 1;
    varargout(countOut) = {x_estEKS}; countOut = countOut + 1;
    varargout(countOut) = {x_errVarEKS}; countOut = countOut + 1;
    
    countTrack = countTrack + 1;     % Increment counter
end

varargout(countOut) = {x}; countOut = countOut + 1;

if plot_flag
    %% Initial Plotting Variables
    titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
    titleCell(2,1)  = {'r'};            % Color for true state plot

    % A basic plotting routine to visualize results
    for j = 1:sum(algFlag)
        plotStateTracksFZ(trueState,estTracks(:,:,j),titleCell(:,[1 j+1]), nP);
    end
end