function [] = checkTriangle()

freqs = [550 900 1500];
bws   = [100 100 100];

fs   = 10000; % Sampling rate

arCoeffs = [1];
for i = 1:length(freqs)
    r = exp(-pi*bws(i)/fs); % Compute BW
    % Get AR Coefficients
    a(1) = -2*r*cos(pi*freqs(i)/(fs/2));
    a(2) = r^2;
    arCoeffs = conv(arCoeffs, [1 a]);
end

figure; freqz(1, arCoeffs,512,fs);

cepOrder = 10;
C  = lpc2c(-arCoeffs(2:end)',cepOrder);
C2 = fb2cp(freqs,bws,cepOrder,fs);
disp([C C2']);
dd = 3;