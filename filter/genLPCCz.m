function C = genLPCCz(wav, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder)

% Function generates LPC cepstral coefficients from waveform data,
% estimating poles and zeros
% 
% INPUT
%   wav - Waveform
%   win - Window to use
%   wOverlap - In percent, how much to overlap (determines frame shift)
%   peCoeff  - Pre-emphasis coeffient
%   lpcOrder - Number of LPC Coefficients to use
%   zOrder - Number of MA Coefficients to use (two times anti-resonances)
%   cepOrder - Number of Cepstral Coefficients to use
%
% OUTPUT
%   C - LPCC coeffients (observations)
%
% Author: Daryush
% Created:  4/23/10
% Modified: 4/23/10

% Pre-emphasize using a first order difference (1-zero filter)
wav = filter([1 -peCoeff],1,wav);

wLength = length(win); % Compute window length (in samples)

% Truncate samples that don't fall within a fixed number of windows
% (PW: Seems this is what VTR database/Wavesurfer does)
wav = wav(1:end-mod(length(wav),wLength)); %DR: this might be an issue
sLength = length(wav);

% Compute number of frames
numFrames = floor(sLength/((1-wOverlap)*wLength))-1;

% Compute left and right boundaries
wLeft  = 1:wLength*(1-wOverlap):sLength-wLength+1;
wRight = wLength:wLength*(1-wOverlap):sLength;

% Store all extracted LPC coefficients
allCoeffsP = zeros(numFrames, lpcOrder+1);
allCoeffsZ = zeros(numFrames, zOrder+1);

for i=1:numFrames
    % Pull out current segment and multiply it by window
    curSegment = wav(wLeft(i):wRight(i));
    
    % Estimate ARMA parameters using armax function from Sys. ID. toolbox
    data = iddata(win.*curSegment,[],1); % Package input
    m = armax(data,[lpcOrder zOrder]); % Call estimator with desired model orders
    
    allCoeffsP(i,:) = m.a;
    allCoeffsZ(i,:) = m.c;
end

% Convert ARMA coefficients to cepstral coefficients
C1 = lpc2c(-allCoeffsP(:,2:end)',cepOrder);
C2 = lpc2c(-allCoeffsZ(:,2:end)',cepOrder);

C = C1 - C2;