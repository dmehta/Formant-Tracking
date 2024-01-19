function C = genCeps(wav, win, wOverlap, peCoeff, cepOrder)

% Function generates cepstral coefficients from real or complex cepstrum
% directly from waveform data.
% 
% INPUT
%   wav - Waveform
%   win - Window to use
%   wOverlap - In percent, how much to overlap (determines frame shift)
%   peCoeff  - Pre-emphasis coeffient
%   cepOrder - Number of Cepstral Coefficients to use
%
% OUTPUT
%   C - coefficients from cepstrum (observations)

% Author: Daryush
% Created:  11/28/10
% Modified: 11/28/10

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

% Store first N coefficients of complex cepstrum
C = zeros(cepOrder, numFrames);

for i=1:numFrames
    % Pull out current segment and multiply it by window
    curSegment = wav(wLeft(i):wRight(i));
    Ctemp = rceps(win.*curSegment); % rceps or cceps
    C(:,i) = Ctemp(2:cepOrder+1,1); % first coef is 0th so start with second
end