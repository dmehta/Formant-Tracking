clear
close all

addpath(genpath('../')); % Set paths

fs = 16e3; % Hz
cepOrder = 15;

%% set center frequencies and bandwidths
% F = []; Fbw = [];
% F = [500 1500 2500 3500]; Fbw = 80+40*[0:3];
% F = [900 1500]; Fbw = 80+40*[0:1];
F = 500; Fbw = 100;
F=F'; Fbw = Fbw';

% Z = []; Zbw = [];
Z = 2000; Zbw = 80;
% Z = [500 1500 2500 3500]+100; Zbw = 80+40*[0:3];
Z=Z'; Zbw = Zbw';

%% route 1 vs route 2
% [num, denom] = fb2tf(F, Fbw, Z, Zbw, fs);
% figure, freqz(num, denom, 512, fs)
% C = lpc2cz(-denom(2:end)',-num(2:end)',cepOrder); % for speech to LPCC coefficients
% 
% % route 2
% C2 = fb2cpz(F, Fbw, Z, Zbw, cepOrder, fs); % for estimation equation
% 
% [C C2']
% sum(C-C2')

%% Create an ARMA model by filtering a white noise sequence
dur = 0.25; % in s
N = round(dur*fs);
[num, denom] = fb2tf(F, Fbw, Z, Zbw, fs);
x = filter(num, denom, randn(N,1));

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
legend('Transfer function', 'Periodogram')

%%
% Estimate AR params using arcov (this leads to a biased estimate)
[arCoeffs e] = arcov(x,length(F)*2);
disp(' ')
disp('True Coefficients');
disp(['AR Coeffs:' num2str(denom)])
disp(['MA Coeffs:' num2str(num)])
disp('Covariance Method: Estimated AR Coefficients')
disp(['AR Coeffs:' num2str(arCoeffs)])
[spec, freq] = freqz(1, arCoeffs, 512, fs);
figure, subplot(211)
plot(freq, 20*log10(abs(spec)), 'b')
[spec, freq] = freqz(1, denom, 512, fs);
hold on, subplot(211)
plot(freq, 20*log10(abs(spec)), 'r')
title('AR estimate spectrum (ARCOV)')
legend('Estimate', 'True')

%% Estimate ARMA parameters using armax function from Sys. ID. toolbox
data = iddata(x,[],1); % Package input
m = armax(data,[length(F)*2 length(Z)*2]); % Call estimator with desired model orders
m.a = denom;
m.c = num;
disp('Sys ID toolbox ARMA estimates');
disp(['AR Coeffs: ' num2str(m.a)]); % Estimated AR Coefficients
disp(['MA Coeffs: ' num2str(m.c)]); % Estimated MA Coefficients

[spec, freq] = freqz(1, m.a, 512, fs);
figure, subplot(211)
plot(freq, 20*log10(abs(spec)), 'b')
[spec, freq] = freqz(1, denom, 512, fs);
hold on, subplot(211)
plot(freq, 20*log10(abs(spec)), 'r')
title('AR estimate spectrum (ARMA)')
legend('Estimate', 'True')

[spec, freq] = freqz(1, m.a, 512, fs);
figure, subplot(211)
plot(freq, 20*log10(abs(spec)), 'b')
[spec, freq] = freqz(1, denom, 512, fs);
hold on, subplot(211)
plot(freq, 20*log10(abs(spec)), 'r')
title('AR estimate spectrum (ARMA)')
legend('Estimate', 'True')

[spec, freq] = freqz(m.c, m.a, 512, fs);
figure, subplot(211)
plot(freq, 20*log10(abs(spec)), 'b')
[spec, freq] = freqz(num, denom, 512, fs);
hold on, subplot(211)
plot(freq, 20*log10(abs(spec)), 'r')
title('ARMA estimate spectrum (ARMA)')
legend('Estimate', 'True')

%% Calculate LPCC coefficients in each window
wType = 'hamming';  % window type
wLengthMS  = 20;    % Length of window (in milliseconds)
wOverlap = 0.5;     % Factor of overlap of window
lpcOrder = length(F)*2;       % Number of LPC coefficients
zOrder = length(Z)*2;         % Number of MA coefficients
peCoeff = 0;     % Pre-emphasis factor

wLength = floor(wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(wType,wLength);

y = genLPCCz(x, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder);

%% Do a plot of the observations
% figure;
% imagesc(log(abs(y))); colorbar;
% title('Cepstral Coefficients');
% xlabel('Frame Number');

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

SNR = 25; % dB, can't be too high
[y, oNoiseVar] = addONoise(y, SNR);
R = oNoiseVar*eye(cepOrder);       % Measurement noise covariance matrix R

bwFlag = 1; % 0 - Use loaded bandwidths, 1 - Average bandwidths
bwStates = repmat([Fbw; Zbw], 1, size(y,2));

% A voice activity detector is not used here in the synthetic case
formantInds = ones(N,nP + nZ);

% General Settings
algFlag = [1 0]; % Select 1 to run, 0 not to
EKF = 1; EKS = 2;

% Initial state of formant trackers
% x0 = trueState(:,1);
x0 = [F; Z]+500;

countTrack = 1; % Counter for storing results
% Initialize root-mean-square error matrices:
rmse    = zeros(nP + nZ, sum(algFlag));
relRmse = zeros(nP + nZ, sum(algFlag));

% Run Extended Kalman Filter
if algFlag(EKF)
    smooth = 0;

    [x_estEKF x_errVarEKF] = formantTrackEKSZ(y, Fmatrix, Q, R, x0, formantInds, fs, bwStates, nP, smooth);
    %[x_estEKF x_errVarEKF] = formantTrackEKSZ(y, F, Q, R, x0, formantInds, fs, bwStates, nP, smooth);

    %Track estimate into data cube for plot routines
    estTracks(:,:,countTrack) = x_estEKF;
    estVar(:,:,:,countTrack) = x_errVarEKF;
    titleCell(1,countTrack+1) = {'EKF'};
    titleCell(2,countTrack+1) = {'g-.'};

    countTrack = countTrack + 1;     % Increment counter
    figure, plot(x_estEKF(1,:))
    title('EKF')
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
    title('EKS')
end

%% Initial Plotting Variables
break
titleCell(1,1)  = {'True State'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};            % Color for true state plot

% A basic plotting routine to visualize results
plotStateTracksFZ(trueState,estTracks(:,:,1),titleCell(:,[1 2]), nP);
plotStateTracksFZ(trueState,estTracks(:,:,2),titleCell(:,[1 3]), nP);