clear
close all

addpath(genpath('../')); % Set paths

fs = 10e3; % Hz
cepOrder = 15;

%% set center frequencies and bandwidths
% F = []; Fbw = [];
% F = [500 1500 2500 3500]; Fbw = 80+40*[0:3];
F = [900 1500]; Fbw = 80+40*[0:1];
F=F'; Fbw = Fbw';

% Z = []; Zbw = [];
% Z = 950; Zbw = 80;
Z = [500 1500 2500 3500]+100; Zbw = 80+40*[0:3];
Z=Z'; Zbw = Zbw';

%% Create an ARMA model by filtering a white noise sequence

% Compute transfer function coefficients
[num, denom] = fb2tf(F, Fbw, Z, Zbw, fs);
[spec, freq] = freqz(num, denom, 512, fs);
figure, freqz(num, denom, 512, fs)

dur = 0.25; % in s
N = round(dur*fs);
x = filter(num, denom, randn(N,1));

figure, subplot(211)
plot(x)
title('Time-domain waveform');
xlabel('Samples'); ylabel('Amplitude');

hold on, subplot(212)
[spec_welch, freq_welch] = pwelch(x-mean(x), [], [], [], fs);
spec_welch = 10*log10(spec_welch);
plot(freq_welch, spec_welch)
title('Power spectrum');
xlabel('Frequency (Hz)'); ylabel('Power (dB)');

hold on, subplot(212)
plot(freq, 20*log10(abs(spec))-38, 'r')
legend('Waveform', 'True')

%%
% Estimate AR params using arcov (this leads to a biased estimate)
[arCoeffs e] = arcov(x,length(F)*2);
disp(' ')
disp('True Coefficients');
disp(['AR Coeffs:' num2str(denom/denom(1))])
disp(['MA Coeffs:' num2str(num/num(1))])
disp('Covariance Method: Estimated AR Coefficients')
disp(['AR Coeffs:' num2str(arCoeffs)])
[spec, freq] = freqz(1, arCoeffs, 512, fs);
figure, subplot(211)
plot(freq, 20*log10(abs(spec)), 'b')

[spec, freq] = freqz(1, denom, 512, fs);
hold on, subplot(211)
plot(freq, 20*log10(abs(spec)), 'r')

%% Estimate ARMA parameters using armax function from Sys. ID. toolbox
data = iddata(x,[],1); % Package input
m = armax(data,[length(F)*2 length(Z)*2]); % Call estimator with desired model orders
disp('Sys ID toolbox ARMA estimates');
disp(['AR Coeffs: ' num2str(m.a)]); % Estimated AR Coefficients
disp(['MA Coeffs: ' num2str(m.c)]); % Estimated MA Coefficients
[spec, freq] = freqz(1, m.a, 512, fs);
figure, subplot(211)
plot(freq, 20*log10(abs(spec)), 'b')

[spec, freq] = freqz(1, denom, 512, fs);
hold on, subplot(211)
plot(freq, 20*log10(abs(spec)), 'r')

[spec, freq] = freqz(m.c, m.a, 512, fs);
figure, subplot(211)
plot(freq, 20*log10(abs(spec)), 'b')

[spec, freq] = freqz(num, denom, 512, fs);
hold on, subplot(211)
plot(freq, 20*log10(abs(spec)), 'r')

%% Calculate LPCC coefficients in each window
wType = 'hamming';  % window type
wLengthMS  = 20;    % Length of window (in milliseconds)
wOverlap = 0.5;     % Factor of overlap of window
lpcOrder = 4;       % Number of LPC coefficients
zOrder = 8;         % Number of MA coefficients
peCoeff = 0;     % Pre-emphasis factor

wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);

y = genLPCCz(x, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);

%% Do a plot of the observations
figure;
imagesc(log(abs(y))); colorbar;
title('Cepstral Coefficients');
xlabel('Frame Number');

%% %%%%%%%%%%%%%%%%%%% Tracking Algorithms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for tracking algorithms %%%
% x_{k+1} = Fx_{k} + w_k, w_k ~ N(0, Q)
% y_{k}   = Hx_{k} + v_k, v_k ~ N(0, R)
% We need to set the parameters: F, Q and R
% H is obtained in the EKF via linearization about the state

nP = length(F);
nZ = length(Z);

Fmatrix = eye(nP + nZ);   % Process Matrix F

qFormantVar = 100;
Q = diag(ones(nP+nZ,1)*qFormantVar);

SNR = 20; % dB
[y, oNoiseVar] = addONoise(y, SNR);
R = oNoiseVar*eye(cepOrder);       % Measurement noise covariance matrix R

bwFlag = 1; % 0 - Use loaded bandwidths, 1 - Average bandwidths
bwStates = repmat([Fbw; Zbw], 1, size(y,2));

% A voice activity detector is not used here in the synthetic case
formantInds = ones(N,nP + nZ);

% General Settings
algFlag = [0 1]; % Select 1 to run, 0 not to
EKF = 1; EKS = 2;

% Initial state of formant trackers
% x0 = trueState(:,1);
x0 = [F; Z];

countTrack = 1; % Counter for storing results
% Initialize root-mean-square error matrices:
rmse    = zeros(nP + nZ, sum(algFlag));
relRmse = zeros(nP + nZ, sum(algFlag));

% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;
    [x_estEKF x_errVarEKF] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);

    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'g-.'};

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

    figure, plot(x_estEKS')
end
break
%% Initial Plotting Variables
titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};            % Color for true state plot

% A basic plotting routine to visualize results
plotStateTracksFZ(trueState,estTracks(:,:,1),titleCell(:,[1 2]), nP);
plotStateTracksFZ(trueState,estTracks(:,:,2),titleCell(:,[1 3]), nP);