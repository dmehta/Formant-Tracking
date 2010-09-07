addpath(genpath('../')); % Paths

%% parameters

% synthesis parameters
N = 100; % number of observations
fs = 8e3; % in Hz
rand('state',2); randn('state',2); % Seeds for JASA figure; and delete any rand('state') in runSynthZ

F = [500]'; % in Hz, time-invariant formant center frequency (F)
Fbw = [50]'; % in Hz, time-invariant formant bandwidth (Fbw)
Z = [700]'; % in Hz, time-invariant formant center frequency (Z)
Zbw = [50]'; % in Hz, time-invariant formant bandwidth (Zbw)

pNoiseVar(1) = 1000;
snr_dB(1) = 40;
cepOrder(1) = 30; %max(aParams.lpcOrder, aParams.zOrder);

% tracker parameters
pNoiseVar(2) = 1000;
snr_dB(2) = 20;
cepOrder(2) = 15; %max(aParams.lpcOrder, aParams.zOrder);
trackBW = 1;
plot_flag = 0; % do plots
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
offset = 100; % set initial state offset, in Hz

%% misc calculations
if trackBW
    x0 = [F; Fbw; Z; Zbw]+offset;
else
    x0 = [F; Z]+offset;
end

%% synthesize and track
[rmse, x_est, x_errVar, trueState] = runSynthZ(F, Fbw, Z, Zbw, N, pNoiseVar, snr_dB, ...
        cepOrder, fs, trackBW, plot_flag, algFlag, x0);
    
%% plot tracks individually in grid style with covariances
titleCell(1,2) = {'EKS'}; % hard coded for now
titleCell(2,2) = {'b:'};
titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};      % Color for true state plot
nP = size(F,1);

plotStateTracksFZ_EstVar_Truth(trueState,x_est,x_errVar,titleCell,nP,trackBW)
disp(['Mean RMSE: ', num2str(mean(rmse))])
rmse