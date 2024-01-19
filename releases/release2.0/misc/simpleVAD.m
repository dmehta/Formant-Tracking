% Simple voice activity detection
%
% INPUT
%
% analysisParams - Datastructure specifying the short-time analysis scheme
%                  within which to do detection
%
% quantThresh    - Simple quantile threshold on power
% multiBand      - Boolean indicating whether to do VAD in more than one band
% doPlot         - Show the mask that is produced
% Optional Params -
%
%
% Author: Daniel Rudoy
% Created:  12/15/2007
% Modified: 12/16/2007
function [frameInds] = simpleVAD(wav, analysisParams, quantThresh, multiBand, bandEdges, doPlot)

% Extract the necessary parameters from the analysis set
fs       = analysisParams.fs;         % Sampling rate
wLength  = analysisParams.wLength;    % Window length
wOverlap = analysisParams.wOverlap;   % Window overlap
win      = analysisParams.win;        % Actual window
peCoeff  = analysisParams.peCoeff;    % Pre-emphasis coefficient
wLengthMS = analysisParams.wLengthMS;

% Pre-emphasize waveform
wav = filter([1 -peCoeff],[1],wav);

if(~multiBand)
    display('Single-band voice activity detection');
    bandEdges = [1 analysisParams.fs/2];
else
    display('Multi-band voice activity detection');

    % Make sure does not exceed allotted bandwidth
    bandEdges((bandEdges > fs/2)) = fs/2;
end

numBands = size(bandEdges,1);

% Truncate samples that don't fall within a fixed number of windows
% (PJW: Seems this is what VTR database/Wavesurfer does)
%wav = wav(1:end-mod(length(wav),wLength)); %DGR: this might be an issue
numFrames = floor(length(wav)/(wLengthMS/1000*fs/2))-1;
sLength = round((numFrames+1)*wLengthMS/1000*fs/2);
wav = wav(1:sLength);

% Compute number of frames
% numFrames = floor(sLength/((1-wOverlap)*wLength))-1;

% Compute left and right boundaries
% This is only checked for wOverlap = 1/2, boundary problems may exist
wLeft  = 1:wLength*(1-wOverlap):sLength-wLength+1;
wRight = wLength:wLength*(1-wOverlap):sLength;

% Pre-allocate array
energy = zeros(numFrames, numBands);

for i=1:numFrames
    % Pull out current segment and multiply it by window
    curSegment = wav(wLeft(i):wRight(i));
    winSegment = win.*curSegment;

    % Compute energy for each energy band
    psd = abs(fft(winSegment));
    psdRight = psd(1:wLength/2).^2;
    energyLocsInd = ceil(bandEdges./(fs/2).*length(psdRight));
    for j = 1:numBands
        energy(i,j) = sum(psdRight(energyLocsInd(j,1):energyLocsInd(j,2)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do thresholding in each band (if there is one band, just do it across the
% whole thing)
frameInds = ones(numBands,numFrames);
for j = 1:numFrames
    inds(j,:) = energy(j,:) < quantile(energy,[quantThresh]); % #ok<AGROW>
end
% Zero out the indices where coasting should occur
frameInds(inds) = 0;

if(doPlot)
    figure;
    imagesc(frameInds);
    title('Result of Simple Voice Activity Detector');
    xlabel('Short-time frame number');
end







