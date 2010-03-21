function [] = testARMAX()

% Play around with the SysID toolbox to see about pole zero estimation
close all; clear all;

% Synthesize ARMA model

N = 1000; % Number of data points
p = 2; % Number of poles 
q = 2; % Number of zeros


rz = .9; thetaz =  pi/4; % Complex-conjugate zero pair:  (+/-  pi/4, .9)
rp = .99; thetap = pi/3;  % Complex-conjugate pole pair: (+/- pi/3, .99)

% Create an ARMA model by filtering a white noise sequence
x = filter([1 -2*cos(thetaz)*rz rz^2], [1 -2*cos(thetap)*rp rp^2], randn(N,1));

disp(['True AR Coefficients:' num2str([1, -2*cos(thetap)*rp, rp^2])]);
disp(['True MA Coefficients:' num2str([1, -2*cos(thetaz)*rz, rz^2])]);

figure;
subplot(2,1,1);
plot(x);
title('Time-domain waveform');
xlabel('Samples'); ylabel('Amplitude');
subplot(2,1,2);

lfft = 1024;
X = 20*log10(abs(fft(x,lfft)));
X = X(1:1:lfft/2);
bins = linspace(0,pi,length(X));
plot(bins, X);
title('Log-Power Spectrum');
xlabel('Frequency (Rad)'); ylabel('Magnitude (dB)');
axis('tight');

% Estimate AR params using arcov (this leads to a biased estimate)
[arCoeffs e] = arcov(x,p);
disp('Covariance Method: Estimated AR Coefficients');
disp(['AR Coeffs:' num2str(arCoeffs)]);

% Estimate ARMA parameters using armax function from Sys. ID. toolbox
data = iddata(x,[],1); % Package input
m = armax(data,[p q]); % Call estimator with desired model orders
disp('Sys ID toolbox ARMA estimates');
disp(['AR Coeffs: ' num2str(m.a)]); % Estimated AR Coefficients
disp(['MA Coeffs: ' num2str(m.c)]); % Estimated MA Coefficients
