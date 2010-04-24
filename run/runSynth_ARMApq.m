%% runSynth_ARMApq
clear

fs = 10e3; % Hz
cepOrder = 15;
dur = 0.25; % in s
N = round(dur*fs);

% set center frequencies and bandwidths
% F = []; Fbw = [];
% F = [500 1500 2500 3500]; Fbw = 80+40*[0:3];
F = [900 1500]; Fbw = 80+40*[0:13];
F=F'; Fbw = Fbw';

% Z = []; Zbw = [];
Z = 1000; Zbw = 100;
% Z = [500 1500 2500 3500]+100; Zbw = 80+40*[0:3];
% Z=Z'; Zbw = Zbw';

[num, denom] = fb2tf(F, Fbw, Z, Zbw, fs);
[spec, freq] = freqz(num, denom, 512, fs);
% figure, freqz(num, denom, 512, fs)


%% Create an ARMA model by filtering a white noise sequence
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

% Estimate ARMA parameters using armax function from Sys. ID. toolbox
data = iddata(x,[],1); % Package input
m = armax(data,[length(F)*2 length(Z)*2]); % Call estimator with desired model orders
disp('Sys ID toolbox ARMA estimates');
disp(['AR Coeffs: ' num2str(m.a)]); % Estimated AR Coefficients
disp(['MA Coeffs: ' num2str(m.c)]); % Estimated MA Coefficients

%%
% % C = lpc2cz(-denom(2:end)',-num(2:end)',cepOrder) % for speech to LPCC coefficients
% C = lpc2c(-denom(2:end),cepOrder) % for speech to LPCC coefficients
% 
% %%
% % C2 = fb2cpz(F, Fbw, Z, Zbw, cepOrder, fs); % for estimation equation
% % C2'
% 
% C2 = fb2cp(F, Fbw, cepOrder, fs); % for estimation equation
% C2'
% 
% %%

