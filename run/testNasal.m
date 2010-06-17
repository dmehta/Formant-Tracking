%% parameters

% WAV file to analyze
% filename = '../data/DDM_speech/WAV/m.wav';
filename = '../data/DDM_speech/WAV/n.wav';
% filename = '../data/DDM_speech/WAV/ng.wav';
% filename = '../data/DDM_speech/WAV/m_asp.wav';
% filename = '../data/DDM_speech/WAV/n_asp.wav';
% filename = '../data/DDM_speech/WAV/ng_asp.wav';
% filename = '../data/DDM_speech/WAV/an.wav';

% analysis parameters
fs = 10000; % sampling rate (in Hz) to resample to
pOrder = 24;
zOrder = 4;
peCoeff = .9;
dur = 50e-3; % window duration, in s

%%
[x, fs_in] = wavread(filename);
x = resample2fs(x, fs, fs_in);
x = filter([1 -peCoeff],1,x);

win = hamming(round(dur*fs))';
x = win.*x(1:round(dur*fs));

%%
data = iddata(x',[],1); % Package input
m = armax(data, [pOrder zOrder]);

% Estimate AR params using arcov (this leads to a biased estimate)
[arCoeffs e] = arcov(x,pOrder);

%%
fftlen = 2^(nextpow2(dur*fs));
figure, hold on

subplot(211)
plot(x)
title('Time-domain waveform');
xlabel('Samples'); ylabel('Amplitude');

subplot(212), hold on
spec = 20*log10(abs(fft(x, fftlen)));
freq = (0:length(spec)-1)/length(spec)*fs;
plot(freq(1:end/2), spec(1:end/2), 'Color', [.8 .8 .8])

[spec, freq] = freqz(1, arCoeffs, fftlen, fs);
plot(freq, 20*log10(abs(spec))-25, 'k', 'LineWidth', 2)

[spec, freq] = freqz(m.c, m.a, fftlen, fs);
plot(freq, 20*log10(abs(spec))-25, 'r', 'LineWidth', 2)

ylim([-60 20])
title('Power spectrum');
xlabel('Frequency (Hz)'); ylabel('Power(dB)');
legend('Data spectrum', 'AR', 'ARMA')

%% frequencies of zeros and poles
[b, a, nn, mm] = eqtflength(m.c, m.a);
[z, p, k] = tf2zp(b, a);

ii = find(imag(z) == 0);
z(ii) = [];
z = z(1:2:end);
ii = find(imag(p) == 0);
p(ii) = [];
p = p(1:2:end);

F = angle(p)/(2*pi)*fs;
FBW = -log(abs(p))*fs/pi; % needs to be CHANGED IN FORMANTS.PDF, add negative sign in (1)
ii = find(FBW < 300);
F = F(ii);
FBW = FBW(ii);
[F, jj] = sort(F);
FBW = FBW(jj);

antiF = angle(z)/(2*pi)*fs;
antiFBW = -log(abs(z))*fs/pi; % needs to be CHANGED IN FORMANTS.PDF, add negative sign in (1)
ii = find(antiFBW < 300);
antiF = antiF(ii);
antiFBW = antiFBW(ii);
[antiF, jj] = sort(antiF);
antiFBW = antiFBW(jj);

disp(['Resonance: ', num2str(F')])
disp(['Resonance BW: ', num2str(FBW')])
disp(['Anti-Resonance: ', num2str(antiF')])
disp(['Anti-Resonance BW: ', num2str(antiFBW')])

%% residual
e_arma = filter(m.a, m.c, x);
e_arma_var = mean(e_arma.^2);
disp(['Variance of ARMA residual: ', num2str(e_arma_var)])

e_ar = filter(arCoeffs, 1, x);
e_ar_var = mean(e_ar.^2);
disp(['Variance of AR residual: ', num2str(e_ar_var)])