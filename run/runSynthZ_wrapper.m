%% parameters

% synthesis parameters
N = 100; % number of observations
fs = 8e3; % in Hz

F = [500 1000]'; % in Hz, time-invariant formant center frequency (F)
Fbw = [50 50]'; % in Hz, time-invariant formant bandwidth (Fbw)
Z = [700 1400]'; % in Hz, time-invariant formant center frequency (Z)
Zbw = [100 50]'; % in Hz, time-invariant formant bandwidth (Zbw)

% analysis parameters
aParams.lpcOrder = length(F)*2; % Number of AR coefficients
aParams.zOrder = length(Z)*2;   % Number of MA coefficients
aParams.fs = fs;                % sampling rate (in Hz)

% tracker parameters
pNoiseVar = 5;
snr_dB = 100;
cepOrder = 30; %max(aParams.lpcOrder, aParams.zOrder);
trackBW = 1;
plot_flag = 0; % do plots
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
offset = 0; % set initial state offset, in Hz

%% misc calculations
if trackBW
    x0 = [F; Fbw; Z; Zbw]+offset;
else
    x0 = [F; Z]+offset;
end

%% synthesize and track
[rmse, x_est, x_errVar, trueState] = runSynthZ(F, Fbw, Z, Zbw, N, pNoiseVar, snr_dB, ...
        cepOrder, fs, trackBW, plot_flag, algFlag, x0, aParams);
    
%% plot tracks individually in grid style with confidence intervals
titleCell(1,2) = {'EKS'}; % hard coded for now
titleCell(2,2) = {'b:'};
titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};      % Color for true state plot
nP = size(F,1);

plotStateTracksFZ(trueState,x_est,titleCell(:,[1 2]), nP, trackBW);
disp(['Mean RMSE: ', num2str(mean(rmse))])