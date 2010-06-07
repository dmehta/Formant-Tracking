%%
clear 

fs = 10000; % Hz
pOrder = 12;
zOrder = 12;
peCoeff = .9;
dur = 40e-3; % s

filename = '../data/DDM_speech/WAV/n.wav';

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
figure, hold on

subplot(211)
plot(x)
title('Time-domain waveform');
xlabel('Samples'); ylabel('Amplitude');

subplot(212), hold on

spec = 20*log10(abs(fft(x, 512)));
freq = (0:length(spec)-1)/length(spec)*fs;
plot(freq(1:end/2), spec(1:end/2), 'Color', [.8 .8 .8])

[spec, freq] = freqz(sqrt(e), arCoeffs, 512, fs);
plot(freq, 20*log10(abs(spec))+20, 'k', 'LineWidth', 2)

[spec, freq] = freqz(m.c, m.a, 512, fs);
plot(freq, 20*log10(abs(spec))-20, 'r', 'LineWidth', 2)

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
e_arma_var = mean(e_arma.^2)

e_ar = filter(arCoeffs, 1, x);
e_ar_var = mean(e_ar.^2)

%%
%                 % focus on changing first three formants, so find the three poles
%                 % closest to the unit circle, not counting real poles
%                 % corresponding to DC (freq = 0)
%                 i = find(imag(z) == 0);
%                 newZeros = z(i);
% 
%                 z(i) = -2;
% 
%                 [y, i] = (max(real(z)));
%                 Z1 = z(i);
%                 z(i) = -2; z(i+1) = -2; % assuming complex conjugate is next
% 
%                 %                 [y, i] = (max(real(z)));
%                 %                 Z2 = z(i);
%                 %                 z(i) = -2; z(i+1) = -2; % assuming complex conjugate is next
%                 % 
%                 %                 [y, i] = (max(real(z)));
%                 %                 Z3 = z(i);
%                 %                 z(i) = -2; z(i+1) = -2; % assuming complex conjugate is next        
%                 % 
%                 %                 newZeros = [newZeros; z( find( z ~= -2) )]; % keep other zeroes same
% 
%                 antiF = angle(Z1)/(2*pi)*fs;
%                 antiFBW = -log(abs(Z1))*fs/pi; % needs to be CHANGED IN FORMANTS.PDF, add negative sign in (1)
                
%                 % focus on changing first three formants, so find the three poles
%                 % closest to the unit circle, not counting real poles
%                 % corresponding to DC (freq = 0)
%                 i = find(imag(p) == 0);
%                 newPoles = p(i);
% 
%                 p(i) = -2;
% 
%                 [y, i] = (max(real(p)));
%                 F1 = p(i);
%                 p(i) = -2; p(i+1) = -2; % assuming complex conjugate is next
% 
%                 [y, i] = (max(real(p)));
%                 F2 = p(i);
%                 p(i) = -2; p(i+1) = -2; % assuming complex conjugate is next
% 
%                 %                 [y, i] = (max(real(p)));
%                 %                 F3 = p(i);
%                 %                 p(i) = -2; p(i+1) = -2; % assuming complex conjugate is next        
%                 
%                 newPoles = [newPoles; p( find( p ~= -2) )]; % keep other poles same
% 
%                 F = angle([F1 F2])/(2*pi)*fs;
%                 FBW = -log(abs([F1 F2]))*fs/pi; % needs to be CHANGED IN FORMANTS.PDF, add negative sign in (1)