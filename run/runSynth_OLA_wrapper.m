% addpath(genpath('../')); % Paths

%%
clear all

%% parameters
% synthesis parameters
dur = 2; % in s
fs = 10000; % in Hz
sParams.wType = 'hanning';      % window type
sParams.wLengthMS  = 20;       % Length of window (in milliseconds)
sParams.wOverlap = 0.5;         % Factor of overlap of window
snr_dB(1) = 20;
cepOrder(1) = 20; %max(aParams.lpcOrder, aParams.zOrder);
randn('state', 5)

% linear trajectory
Fbeg = [500 2000 2510 3400]'; Fend = [1800 500 2000 3500]';
Fbwbeg = [120 80 160 200]'; Fbwend = [120 80 160 200]';

% Fbeg = [850 1100 2810]'; Fend = [310 2790 3310]';
% Fbwbeg = [80 120 160]'; Fbwend = [80 120 160]';

% Fbeg = [500]'; Fend = [900]';
% Fbwbeg = [100]'; Fbwend = [100]';

% Fbeg = [257 1800]'; Fend = [257 1800]';
% Fbwbeg = [32 50]'; Fbwend = [32 50]';

% /n/
% Fbeg = [257 1491]'; Fend = [257 1491]';
% Fbwbeg = [32 100]'; Fbwend = [32 100]';

% /a/
% Fbeg = [850 1100]'; Fend = [850 1100]';
% Fbwbeg = [80 120]'; Fbwend = [80 120]';

% /n a/
% Fbeg = [257 1491]'; Fend = [850 1100]';
% Fbwbeg = [32 100]'; Fbwend = [80 120]';

% /a n/
% Fbeg = [850 1100]'; Fend = [257 1491]';
% Fbwbeg = [80 120]'; Fbwend = [32 100]';

% Zbeg = [900]'; Zend = [500]';
% Zbwbeg = [100]'; Zbwend = [100]';
% Zbeg = [1400]'; Zend = [1400]';
% Zbwbeg = [100]'; Zbwend = [100]';

% /n/
% Zbeg = [1223]'; Zend = [1223]';
% Zbwbeg = [52]'; Zbwend = [52]';

Zbeg = []'; Zend = []';
Zbwbeg = []'; Zbwend = []';

% analysis parameters
aParams.wType = 'hanning';          % window type
aParams.wLengthMS = sParams.wLengthMS;  % Length of window (in milliseconds)
aParams.wOverlap = sParams.wOverlap;    % Factor of overlap of window
aParams.lpcOrder = length(Fbeg)*2;      % Number of AR coefficients
aParams.zOrder = length(Zbeg)*2;        % Number of MA coefficients
aParams.peCoeff = 0;                    % Pre-emphasis factor (0.9, 0.7, etc.)
aParams.fs = fs;                        % sampling rate (in Hz) to resample to

% tracker parameters
snr_dB(2) = 20;
cepOrder(2) = 20; %max(aParams.lpcOrder, aParams.zOrder);
trackBW = 1;
plot_flag = 0; % plot transfer functions
algFlag = [0 1]; % Select 1 to run, 0 not to; [EKF EKS]
offset = 100; % set initial state offset, in Hz
variableParam = 50;

% Monte Carlo analysis parameters
numTrials = 1;

%% misc calculations
wLength = floor(aParams.wLengthMS/1000*fs);
wLength = wLength + (mod(wLength,2)~=0); % Force even
win = feval(aParams.wType,wLength);
N = floor(dur*aParams.fs);
N = N - mod(N,wLength);
numFrames = floor(N/((1-aParams.wOverlap)*wLength))-1;

if ~isempty(Fbeg) && ~isempty(Zbeg)
    Fsinusoid = 0*repmat(cos(2*pi*8*[0:numFrames-1]/numFrames*dur), size(Fbeg, 1), 1);
    F = interp1([1 numFrames]', [Fbeg Fend]', 1:numFrames)'; if length(Fbeg) == 1, F=F'; end;
    F = F + Fsinusoid + 100*randn(size(F));

    Fbwsinusoid = 0*repmat(cos(2*pi*10*[0:numFrames-1]/numFrames*dur), size(Fbwbeg, 1), 1);
    Fbw = interp1([1 numFrames]', [Fbwbeg Fbwend]', 1:numFrames)'; if length(Fbwbeg) == 1, Fbw=Fbw'; end;
    Fbw = Fbw + Fbwsinusoid + 10*randn(size(Fbw));

    Zsinusoid = 0*repmat(cos(2*pi*7*[0:numFrames-1]/numFrames*dur), size(Zbeg, 1), 1);
    Z = interp1([1 numFrames]', [Zbeg Zend]', 1:numFrames)'; if length(Zbeg) == 1, Z=Z'; end;
    Z = Z + Zsinusoid + 100*randn(size(Z));

    Zbwsinusoid = 0*repmat(cos(2*pi*6*[0:numFrames-1]/numFrames*dur), size(Zbwbeg, 1), 1);
    Zbw = interp1([1 numFrames]', [Zbwbeg Zbwend]', 1:numFrames)'; if length(Zbwbeg) == 1, Zbw=Zbw'; end;
    Zbw = Zbw + Zbwsinusoid + 10*randn(size(Zbw));
    
    if trackBW
        x0 = [F(:, 1); Fbw(:, 1); Z(:, 1); Zbw(:, 1)]+offset;
    else
        x0 = [F(:, 1); Z(:, 1)]+offset;
    end
end

if ~isempty(Fbeg) && isempty(Zbeg)
    Fsinusoid = 0*repmat(cos(2*pi*8*[0:numFrames-1]/numFrames*dur), size(Fbeg, 1), 1);
    F = interp1([1 numFrames]', [Fbeg Fend]', 1:numFrames)'; if length(Fbeg) == 1, F=F'; end;
    F = F + Fsinusoid + 10*randn(size(F));

    Fbwsinusoid = 0*repmat(cos(2*pi*10*[0:numFrames-1]/numFrames*dur), size(Fbwbeg, 1), 1);
    Fbw = interp1([1 numFrames]', [Fbwbeg Fbwend]', 1:numFrames)'; if length(Fbwbeg) == 1, Fbw=Fbw'; end;
    Fbw = Fbw + Fbwsinusoid + 10*randn(size(Fbw));

    Z = [];
    Zbw = [];
    
    if trackBW
        x0 = [F(:, 1); Fbw(:, 1)]+offset;
    else
        x0 = F(:, 1)+offset;
    end
end

if isempty(Fbeg) && ~isempty(Zbeg)
    Zsinusoid = 0*repmat(cos(2*pi*7*[0:numFrames-1]/numFrames*dur), size(Zbeg, 1), 1);
    Z = interp1([1 numFrames]', [Zbeg Zend]', 1:numFrames)'; if length(Zbeg) == 1, Z=Z'; end;
    Z = Z + Zsinusoid + 10*randn(size(Z));

    Zbwsinusoid = 0*repmat(cos(2*pi*6*[0:numFrames-1]/numFrames*dur), size(Zbwbeg, 1), 1);
    Zbw = interp1([1 numFrames]', [Zbwbeg Zbwend]', 1:numFrames)'; if length(Zbwbeg) == 1, Zbw=Zbw'; end;
    Zbw = Zbw + Zbwsinusoid + 10*randn(size(Zbw));

    F = [];
    Fbw = [];
    
    if trackBW
        x0 = [Z(:, 1); Zbw(:, 1)]+offset;
    else
        x0 = Z(:, 1)+offset;
    end
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
    fprintf('.')

    % comment if no variable parameter
    %cepOrder = [cepOrder(1), variableParam(ii)];
    snr_dB = [snr_dB(1), variableParam(ii)];
    
    [rmse{jj,ii}, x_est{jj,ii}, x_errVar{jj,ii}, x{jj,ii}, trueState{jj,ii}] = runSynth_OLA(F, Fbw, Z, Zbw, N, snr_dB, ...
            cepOrder, fs, trackBW, plot_flag, algFlag, x0, aParams, sParams);
    end
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

%% Super-impose over a spectrogram
trial = 1;
variableParamNumber = 1;
nZ = size(Z,1);
aParams.wLength = floor(aParams.wLengthMS/1000*fs);

plotSpecTracks2BW(x{trial,variableParamNumber}, x_est{trial,variableParamNumber}, aParams, nZ, trackBW);
firstline = get(get(gca,'Title'),'String');
title({['Trial #', num2str(trial)]; firstline})
axis tight
format_plot
set(gca, 'PlotBoxAspectRatio', [3 1 1])

%% plot all simulations against truth
figure, hold on
for jj = 1:numTrials
    plot(x_est{jj,variableParamNumber}(1:size(F,1), :)', 'b')
    plot(x_est{jj,variableParamNumber}(size(F,1)+1:end, :)', 'r')
end
plot(trueState{jj,variableParamNumber}', 'LineWidth', 3, 'Color', 'k') % true state is arbitrary!
xlabel('Frame')
ylabel('Frequency (Hz)')
title('Estimated EKS trajectories (Formant-blue, Antiformant-red, True-black)')

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

%% plot tracks individually in grid style with confidence intervals from
%% Monte Carlo analysis - need to fix to pick one trial (repackage x_est)
% titleCell(1,2) = {'EKS'}; % hard coded for now
% titleCell(2,2) = {'b:'};
% titleCell(1,1)  = {'True'};   % Keeps track of trackers used for plotter
% titleCell(2,1)  = {'r'};      % Color for true state plot
% nP = size(F,1);
% 
% plotStateTracksFZ_CI(trueState{trial,variableParamNumber},x_est{trial,variableParamNumber},titleCell(:,[1 2]), nP, trackBW);
% disp(['Mean RMSE: ', num2str(mean(rmse{trial}))])

%% plot frequency vs frame with confidence intervals
% % repackage
% x_estPerFreq = zeros(numTrials, numFrames, size(trueState{jj},1));
% for kk = 1:size(trueState{jj},1)
%     for jj = 1:numTrials
%         x_estPerFreq(jj, :, kk) = x_est{jj}(kk, :);
%     end
% end
% 
% xdata = 1:numFrames;
% figure, box off, hold on
% for kk = 1:size(trueState{jj},1)
%     [ydata_lower ydata_upper ydata] = findCI(x_estPerFreq(:,:,kk), 68);
%     fill([xdata xdata(end:-1:1)], [ydata_lower ydata_upper(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
%     plot(xdata, ydata, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)
%     xlabel('Frame')
%     ylabel('Frequency (Hz)')
% end
% 
% xlim([1 numFrames])
% ylim([0 fs/2])
% format_plot