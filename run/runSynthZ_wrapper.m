% addpath(genpath('../')); % Paths

%% parameters

% synthesis parameters
N = 100; % number of observations
fs = 10e3; % in Hz
% rand('state',sum(100*clock)); randn('state',sum(100*clock)); % Seeds
rand('state',2); randn('state',2); % Seeds for JASA figure

% F = [500 1000]'; % in Hz, time-invariant formant center frequency (F)
% Fbw = [50 50]'; % in Hz, time-invariant formant bandwidth (Fbw)
% Z = [700 1400]'; % in Hz, time-invariant formant center frequency (Z)
% Zbw = [50 50]'; % in Hz, time-invariant formant bandwidth (Zbw)

F = [500 1000]'; % in Hz, time-invariant formant center frequency (F)
Fbw = [50 80]'; % in Hz, time-invariant formant bandwidth (Fbw)
Z = []'; % in Hz, time-invariant formant center frequency (Z)
Zbw = []'; % in Hz, time-invariant formant bandwidth (Zbw)

pNoiseVar(1) = 0;
snr_dB(1) = 50;
cepOrder(1) = 20;

% analysis parameters
cepOrder(2) = 20;

% tracker parameters
pNoiseVar(2) = 0;
snr_dB(2) = 20;
trackBW = 1;
plot_flag = 0; % do plots
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
offset = 0; % set initial state offset, in Hz
variableParam = 50;

% Monte Carlo analysis parameters
numTrials = 1;

%% misc calculations
if trackBW
    x0 = [F; Fbw; Z; Zbw]+offset;
else
    x0 = [F; Z]+offset;
end

%% synthesize and track
rmse = cell(numTrials, 1);
x_est = cell(numTrials, 1);
x_errVar = cell(numTrials, 1);
x = cell(numTrials, 1);
trueState = cell(numTrials, 1);

% loop through
for ii = 1:length(variableParam)
    for jj = 1:numTrials
    %disp(['Processing Trial #', num2str(jj), '...'])
    %fprintf('.')

    % comment if no variable parameter
    %cepOrder = [cepOrder(1), variableParam(ii)];
    snr_dB = [snr_dB(1), variableParam(ii)];
    
    [rmse{jj,ii}, x_est{jj,ii}, x_errVar{jj,ii}, trueState{jj,ii}] = runSynthZ(F, Fbw, Z, Zbw, N, pNoiseVar, ...
            snr_dB, cepOrder, fs, trackBW, plot_flag, algFlag, x0);
    end
    %disp(' ')
end
    
%% plot tracks individually in grid style with covariances from one trial
trial = 1;
variableParamNumber = 1;
titleCell(1,2) = {'EKS'}; % hard coded for now
titleCell(2,2) = {'b:'};
titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
titleCell(2,1)  = {'r'};      % Color for true state plot
nP = size(F,1);

plotStateTracksFZ_EstVar_Truth(trueState{trial,variableParamNumber},x_est{trial,variableParamNumber},x_errVar{trial,variableParamNumber},titleCell,nP,trackBW)
disp(['Mean RMSE: ', num2str(mean(rmse{trial,variableParamNumber}))])
rmse{trial,variableParamNumber}

%% plot RMSE of formant frequency and bandwidth vs variable parameter
xdata = variableParam;

% repackage into 3D matrix of numTrials x length(variableParam) x numStates
rmsePerState = zeros(numTrials, length(variableParam), size(x0,1));
for kk = 1:size(x0,1)
    for jj = 1:numTrials
        for ii = 1:length(variableParam)
            rmsePerState(jj, ii, kk) = rmse{jj,ii}(kk);
        end
    end
end

if trackBW
    rmsePerFreq = zeros(numTrials, length(variableParam), length(F)+length(Z));
    rmsePerBW = rmsePerFreq;
    
    count = 1;
    for kk = [1:length(F), 2*length(F)+1:2*length(F)+length(Z)]
        for jj = 1:numTrials
            for ii = 1:length(variableParam)
                rmsePerFreq(jj, ii, count) = rmse{jj,ii}(kk);
            end
        end
        count = count + 1;
    end
    
    count = 1;
    for kk = [length(F)+1:2*length(F), 2*length(F)+length(Z)+1:2*length(F)+2*length(Z)]
        for jj = 1:numTrials
            for ii = 1:length(variableParam)
                rmsePerBW(jj, ii, count) = rmse{jj,ii}(kk);
            end
        end
        count = count + 1;
    end
else
    rmsePerFreq = rmsePerState;
end
    
% across frequency and bandwidth
stateNumber = 1;
figure, hold on
% [ydata_lower ydata_upper ydata] = findCI(rmsePerState(:,:,stateNumber), 68); % rmse for one state
[ydata_lower ydata_upper ydata] = findCI(mean(rmsePerState, 3), 68); % collapse rmse across all states
fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
% xlabel('# cepstral coefficents')
xlabel('Varied parameter')
ylabel('Average RMSE (Hz)')
format_plot

% across frequency only
figure, hold on
[ydata_lower ydata_upper ydata] = findCI(mean(rmsePerFreq, 3), 68); % collapse rmse across all states
fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
% xlabel('# cepstral coefficents')
xlabel('Varied parameter')
ylabel('Average RMSE of Frequencies (Hz)')
format_plot

% across bandwidth only
figure, hold on
[ydata_lower ydata_upper ydata] = findCI(mean(rmsePerBW, 3), 68); % collapse rmse across all states
fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
% xlabel('# cepstral coefficents')
xlabel('Varied parameter')
ylabel('Average RMSE of Bandwidths (Hz)')
format_plot