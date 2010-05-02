%% make into runSynth_OLA and runAnalysis_OLA when ready

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dur = .5; % in s
snr_dB = 25;
cepOrder = 25;
fs = 16e3;
plot_flag = 1;
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
wType = 'hanning';  % window type
wLengthMS  = 20;    % Length of window (in milliseconds)
wOverlap = 0.5;     % Factor of overlap of window

wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);
N = floor(dur*fs);
N = N - mod(N,wLength);
numFrames = floor(N/((1-wOverlap)*wLength))-1;

%% Piecewise linear trajectory
Fbeg = [500 1500 2500]';
Fend = [1000 1100 3000]';
Fcontour = interp1( [1 numFrames]', [Fbeg, Fend]', 1:numFrames );
Fbw = [100 100 100]';

Zbeg = [4000]';
Zend = [4000]';
Zcontour = interp1( [1 numFrames]', [Zbeg, Zend]', 1:numFrames )';
Zbw = [100]';

% Zcontour = []'; Zbw = []';

%% sinewave modulated trajectory
% Fbeg = [500 1500 2500]';
% Fbw = [100 100 100]';
% Fpert = 75*repmat(cos(2*pi*2*[1:numFrames]/numFrames*dur)', 1, size(Fbeg, 1));
% Fcontour = interp1( [1 numFrames]', [Fbeg Fbeg]', 1:numFrames );
% Fcontour = Fcontour + Fpert;
% 
% Zbeg = [800]';
% Zbw = [100]';
% Zpert = 100*repmat(cos(2*pi*4*[1:numFrames]/numFrames*dur)', 1, size(Zbeg, 1));
% Zcontour = interp1( [1 numFrames]', [Zbeg Zbeg]', 1:numFrames );
% Zcontour = Zcontour' + Zpert;

% Zcontour = []'; Zbw = []';

%% Synthesize each window for overlap-add
x = zeros(N, 1);
windows = zeros(N, 1);

% Compute left and right boundaries
wLeft  = 1:wLength*(1-wOverlap):N-wLength+1;
wRight = wLength:wLength*(1-wOverlap):N;

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
end

x = x./(windows+eps); % divide out window effect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wType = 'hamming'; % window type
wLengthMS  = 20; % Length of window (in milliseconds)
wOverlap = 0.5; % Factor of overlap of window
lpcOrder = size(Fcontour, 2)*2; % Number of LPC coefficients
zOrder = size(Zcontour, 2)*2; % Number of MA coefficients
peCoeff = 0; % Pre-emphasis factor

wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);

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
nP = size(Fcontour, 2);
nZ = size(Zcontour, 2);

Fmatrix = eye(nP + nZ);   % Process Matrix F

pNoiseVar = cov(x(2:end)-x(1:end-1));
Q = diag(ones(nP+nZ,1)*pNoiseVar);

[y, oNoiseVar] = addONoise(y, snr_dB);
R = oNoiseVar*eye(cepOrder); % Measurement noise covariance matrix R

% initial state
if ~isempty(Fcontour) && ~isempty(Zcontour)
    x0 = [Fcontour(1, :)'; Zcontour(1, :)']+100;
end

if ~isempty(Fcontour) && isempty(Zcontour)
    x0 = Fcontour(1, :)'+100;
end

if isempty(Fcontour) && ~isempty(Zcontour)
    x0 = Zcontour(1, :)'+100;
end

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
% trueState = repmat([F; Z], 1, size(y,2));
trueState = [Fcontour'; Zcontour'];

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
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(N);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(N);
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
        rmse(j,countTrack) = norm((estTracks(j,:,countTrack)-trueState(j,:)))/sqrt(N);
        relRmse(j,countTrack) = (rmse(j,countTrack)/norm(trueState(j,:)))*sqrt(N);
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
    titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
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