function [C] = genLPCC(wav, win, wOverlap, peCoeff, lpcOrder, cepOrder);
% LPC Function generation
% INPUT
%   wav - Waveform
%   win - Window to use
%   wOverlap - In percent, how much to overlap (determines frame shift)
%   peCoeff  - Pre-emphasis coeffient
%   lpcOrder - Number of LPC Coefficients to use
%   cepOrder - Number of Cepstral Coefficients to use
%
%
% OUTPUT
%   C - LPCC coeffients (observations)
%
%
% Author: Daniel Rudoy
% Created:  03/13/2007
% Modified: 12/18/2007

% Pre-emphasize using a first order difference (1-zero filter)
wav = filter([1 -peCoeff],1,wav);

% Compute window length (in samples)
wLength = length(win);

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
allCoeffs = zeros(numFrames, lpcOrder+1);

for i=1:numFrames
    % Pull out current segment and multiply it by window
    curSegment = wav(wLeft(i):wRight(i));
    % get the LPC coefficients and associated error
    [allCoeffs(i,:), allError(i)] = lpc(win.*curSegment, lpcOrder);
end

% Convert LPC Coefficients to cepstral coefficients
[C] = lpc2c(-allCoeffs(:,2:end)',cepOrder);
