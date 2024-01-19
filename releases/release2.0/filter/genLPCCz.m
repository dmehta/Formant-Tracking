function C = genLPCCz(wav, win, wOverlap, peCoeff, lpcOrder, zOrder, cepOrder, num, denom)

% Function generates LPC cepstral coefficients from waveform data,
% estimating poles and zeros.
% 
% INPUT
%   wav - Waveform
%   win - Window to use
%   wOverlap - In percent, how much to overlap (determines frame shift)
%   peCoeff  - Pre-emphasis coeffient
%   lpcOrder - Number of AR Coefficients to use
%   zOrder - Number of MA Coefficients to use (two times anti-resonances)
%   cepOrder - Number of Cepstral Coefficients to use
%   num - True coefficients of transfer function numerator (including
%               zeroth coefficient); do not include if not desired
%   denom - True coefficients of transfer function denominator (including
%               zeroth coefficient); do not include if not desired
%
% OUTPUT
%   C - LPCC coeffients (observations)

% Author: Daryush
% Created:  4/23/10
% Modified: 4/30/10 (handle zeros for lpcOrder and zOrder, and truth input)

% Pre-emphasize using a first order difference (1-zero filter)
wav = filter([1 -peCoeff],1,wav);

wLength = length(win); % Compute window length (in samples)

% Truncate samples that don't fall within a fixed number of windows
% (PW: Seems this is what VTR database/Wavesurfer does)
%wav = wav(1:end-mod(length(wav),wLength)); %DR: this might be an issue,
%DDM: it was...now commented out, math down externally. should have divided
%by wLength/2 and subtracted 1 to get number of frames
numFrames = floor(length(wav)/(wLength/2))-1;
sLength = round((numFrames+1)*wLength/2);
wav = wav(1:sLength);

% Compute number of frames
% numFrames = floor(sLength/((1-wOverlap)*wLength))-1;

% Compute left and right boundaries
wLeft  = 1:wLength*(1-wOverlap):sLength-wLength+1;
wRight = wLength:wLength*(1-wOverlap):sLength;

% Store all extracted LPC coefficients
allCoeffsP = zeros(lpcOrder+1, numFrames);
allCoeffsZ = zeros(zOrder+1, numFrames);

for i=1:numFrames
    % Pull out current segment and multiply it by window
    curSegment = wav(wLeft(i):wRight(i));
    
    if exist('num', 'var') && exist('denom', 'var') % ground truth was input
        if iscell(num) && iscell(denom) % frame-by-frame given
            allCoeffsP(:,i) = denom{i}';
            allCoeffsZ(:,i) = num{i}';
        else % one set of coefficients input
            allCoeffsP(:,i) = denom';
            allCoeffsZ(:,i) = num';
        end            
    else
        if zOrder
            % Estimate ARMA parameters using armax function from Sys. ID. toolbox
            data = iddata(win.*curSegment,[],1); % Package input
            m = armax(data,[lpcOrder zOrder]); % Call estimator with desired model orders

            allCoeffsP(:,i) = m.a';
            allCoeffsZ(:,i) = m.c';
        else
            %allCoeffsP(:,i) = arcov(win.*curSegment, lpcOrder)';
            allCoeffsP(:,i) = lpc(win.*curSegment, lpcOrder)';
            allCoeffsZ(:,i) = 1;
        end
    end
end

% Convert ARMA coefficients to cepstral coefficients
C = lpc2cz(-allCoeffsP(2:end,:),-allCoeffsZ(2:end,:),cepOrder);