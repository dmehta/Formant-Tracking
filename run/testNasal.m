%%
clear 

fs = 4000;
pOrder = 4;
zOrder = 4;
peCoeff = .95;

filename = '../data/DDM_speech/WAV/m.wav';

%%
[x, fs_in] = wavread(filename);
x = resample2fs(x, fs, fs_in);
x = filter([1 -peCoeff],1,x);

data = iddata(x',[],1); % Package input
m = armax(data, [pOrder zOrder]);

% Estimate AR params using arcov (this leads to a biased estimate)
[arCoeffs e] = arcov(x,pOrder);

%%
figure, hold on

subplot(211)
plot(-x(1:round(0.05*fs)))
title('Time-domain waveform');
xlabel('Samples'); ylabel('Amplitude');

subplot(212), hold on
ylim([-100 -30])
[spec_true, freq_true] = pwelch(x-mean(x), [], [], [], fs);
spec_true = 10*log10(spec_true);
plot(freq_true, spec_true)

[spec, freq] = freqz(sqrt(e), arCoeffs, 512, fs);
plot(freq, 20*log10(abs(spec))-30, 'k')

[spec, freq] = freqz(m.c, m.a, 512, fs);
plot(freq, 20*log10(abs(spec))-65, 'r')

title('Power spectrum');
xlabel('Frequency (Hz)'); ylabel('Power(dB)');
legend('Peroidogram', 'AR', 'ARMA')

%%
% [b, a, nn, mm] = eqtflength(m.c, m.a);
[z, p, k] = tf2zp(m.c, m.a);

% frequency of zeros
antiF = angle(z)/(2*pi)*fs
antiFBW = -log(abs(z))*fs/pi % needs to be CHANGED IN FORMANTS.PDF, add negative sign in (1)

% frequency of poles
F = angle(p)/(2*pi)*fs
FBW = -log(abs(p))*fs/pi % needs to be CHANGED IN FORMANTS.PDF, add negative sign in (1)


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