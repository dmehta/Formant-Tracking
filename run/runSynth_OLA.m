function varargout = runSynth_OLA(Fcontour, Fbw, Zcontour, Zbw, N, snr_dB, cepOrder, fs, plot_flag, algFlag, x0)

% Track formants and anti-formants (no bandwidths) on synthetic data
% that constructs a waveform using overlap-add of windows generated from an
% ARMA model.
% 
% Author: Daryush Mehta
% Created:  05/09/2010
% Modified: 05/09/2010
% 
% INPUT:
%    Fcontour:  center frequencies of the resonances (numFrames x numFormants), in Hz
%    Fbw:       corresponding bandwidths of the resonances (numFormants x 1), in Hz
%    Z:         center frequencies of the anti-resonances (numFrames x numAntiformants), in Hz
%    Zbw:       corresponding bandwidths of the anti-resonances (numAntiformants x 1), in Hz
%    N:         length of signal, in samples
%    snr_dB:    observation noise, in dB
%    cepOrder:  Number of cepstal coefficients to compute
%    fs:        sampling rate of waveform, in Hz
%    plot_flag: plot figures if 1
%    algFlag:   select 1 to run, 0 not to for [EKF EKS]
%    x0:        initial states (center frequencies) of formant trackers [(numFormants + numAntiformants) x 1], in Hz
% 
% OUTPUT:
%    Depends on algFlag. For each algorithm, two outputs generated--
%       rmse_mean1: average RMSE across all tracks
%       x_est1:  estimated tracks
%       So that if two algorithms run, the following are output:
%       [rmse_mean1, x_est1, rmse_mean2, x_est2]
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% Synthetic Examples:
%   rmse_mean = runSynth_OLA([500; 1000],[100; 100],[],[],.5,200,25,15,16000,1,[1 0],[600; 1600]);
%
% The usual examples from PRAAT and Wavesurfer are not included, because
% it is not clear how to generate paths of zeros along with the paths of
% formants that are available.

addpath(genpath('../')); % Paths
% rand('state',sum(100*clock)); randn('state',sum(100*clock)); % Seeds
rand('state', 2); randn('state', 44);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(N, 1);
windows = zeros(N, 1);
num_all = cell(N, 1);
denom_all = cell(N, 1);

wType = 'hanning';  % window type
wLengthMS  = 20;    % Length of window (in milliseconds)
wOverlap = 0;     % Factor of overlap of window
wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);

% Compute left and right boundaries
wLeft  = 1:wLength*(1-wOverlap):N-wLength+1;
wRight = wLength:wLength*(1-wOverlap):N;

% Compute number of frames
if ~isempty(Fcontour)
    numFrames = size(Fcontour, 1);
else
    numFrames = size(Zcontour, 1);
end

for i=1:numFrames
    if ~isempty(Fcontour) && ~isempty(Zcontour)
        [num, denom] = fb2tf(Fcontour(i, :)', Fbw, Zcontour(i, :)', Zbw, fs);
    end
    
    if ~isempty(Fcontour) && isempty(Zcontour)
        [num, denom] = fb2tf(Fcontour(i, :)', Fbw, [], [], fs);
    end
    
    if isempty(Fcontour) && ~isempty(Zcontour)
        [num, denom] = fb2tf([], [], Zcontour(i, :)', Zbw, fs);
    end
    
    tmp = filter(num, denom, randn(wLength,1));    
    x(wLeft(i):wRight(i)) = x(wLeft(i):wRight(i)) + win.*tmp;
    windows(wLeft(i):wRight(i)) = windows(wLeft(i):wRight(i)) + win;
    
    num_all{i} = num;
    denom_all{i} = denom;
end

x = x./(windows+eps); % divide out window effect

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wType = 'hanning'; % window type
wLengthMS  = 20; % Length of window (in milliseconds)
wOverlap = 0; % Factor of overlap of window
lpcOrder = size(Fcontour, 2)*2; % Number of LPC coefficients
zOrder = size(Zcontour, 2)*2; % Number of MA coefficients
peCoeff = 0; % Pre-emphasis factor

wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);

% y = genLPCCz(x, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder, num_all, denom_all);
y = genLPCCz(x, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);

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
trueState = [Fcontour'; Zcontour'];

nP = size(Fcontour, 2);
nZ = size(Zcontour, 2);
numObs = size(y,2);

Fmatrix = eye(nP + nZ);   % Process Matrix F

pNoiseVar = zeros(1, size(trueState,1));
for ii = 1:length(pNoiseVar)
    pNoiseVar(ii) = cov(trueState(ii, 2:end)-trueState(ii, 1:end-1));
end
Q = diag(pNoiseVar+eps);

[y, oNoiseVar] = addONoise(y, snr_dB);
R = oNoiseVar*eye(cepOrder); % Measurement noise covariance matrix R

bwFlag = 1; % 0 - Use loaded bandwidths, 1 - Average bandwidths
bwStates = repmat([Fbw; Zbw], 1, size(y,2));

% A voice activity detector is not used here in the synthetic case
formantInds = ones(N,nP + nZ);

countTrack = 1; % Counter for storing results
countOut = 1; % Counter for output variables

% Initialize root-mean-square error matrices:
rmse    = zeros(nP + nZ, sum(algFlag));
relRmse = zeros(nP + nZ, sum(algFlag));

EKF = 1; EKS = 2;

%% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;
    [x_estEKF x_errVarEKF] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);

    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'b:'};

    % Compute and Display MSE and RMSE
    for j = 1:nP+nZ
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
    end

    % Display output summary and timing information
    rmse_mean = mean(rmse(:,countTrack));
    varargout(countOut) = {rmse_mean}; countOut = countOut + 1;
    varargout(countOut) = {x_estEKF}; countOut = countOut + 1;
    %display(['Average EKF RMSE: ' num2str(rmse_mean)]);    
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
    for j = 1:nP+nZ
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(numObs);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(numObs);
    end

    % Display output summary and timing information
    rmse_mean = mean(rmse(:,countTrack));
    varargout(countOut) = {rmse_mean}; countOut = countOut + 1;
    varargout(countOut) = {x_estEKS}; countOut = countOut + 1;
    %display(['Average EKS RMSE: ' num2str(rmse_mean)]);
    countTrack = countTrack + 1;     % Increment counter
end

if plot_flag
    %% Initial Plotting Variables
    titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
    titleCell(2,1)  = {'r'};            % Color for true state plot

    % A basic plotting routine to visualize results
    plotStateTracksFZ(trueState,estTracks(:,:,1),titleCell(:,[1 2]), nP);
end

%%
% figure, hold on
% frame = 1:numFrames;
% plot(frame', x_estEKS(1:size(Fcontour, 2), :)', 'b', ...
%     frame', Fcontour(:, 1:size(Fcontour, 2)), 'b:')
% legend('Resonances', 'True resonances')
% xlabel('Frame')
% ylabel('Frequency (Hz)')
% title('Estimated EKS trajectories')

rmse_mean