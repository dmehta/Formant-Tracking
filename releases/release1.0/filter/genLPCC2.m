function[C, formantInds] = genLPCC2(wav, win, wOverlap, peCoeff, lpcOrder, cepOrder, coastJoint, quantThresh, num_formants, fs);
%function[C, formantInds] = genLPCC2(wav, win, wOverlap, peCoeff, lpcOrder, cepOrder, coastJoint, quantThresh, num_formants, fs);
% LPC Function generation
% (with indicator variables)
% INPUT
%   wav - Waveform
%   win - Window to use
%   wType    - Type of window ('hamming', etc)
%   wLength  - Length of window (in samples)
%   wOverlap - In percent, how much to overlap (determines frame shift)
%   peCoeff  - Pre-emphasis coeffient
%   lpcOrder - Number of LPC Coefficients to use
%   cepOrder - Number of Cepstral Coefficients to use
%   coastJoint - Controls whether all coast together or independently
%   quantThresh - Threshold below which energy histogram is truncated
%                 so that a coasting decision can be made
%   num_formants - Number of formants (starting from 1st)
%   fs       - sample rate of argument "wav"

% Formant energy bands (taken from Deng et al, 2006)
fixedLocs = [200 900; 600 2800; 1400 3800; 1700 5000];
energyLocsHz = fixedLocs(1:num_formants,:);
energyLocsHz(find(energyLocsHz>fs/2)) = fs/2;
        
% Pre-emphasize using a first order difference (1-zero filter)
wav = filter([1 -peCoeff],[1],wav);

% Compute window length (in samples)
wLength = length(win);

% Truncate samples that don't fall within a fixed number of windows
% (PW: Seems this is what VTR database/Wavesurfer does)
wav = wav(1:end-mod(length(wav),wLength)); %DR: this might be an issue
sLength = length(wav);

% Compute number of frames
numFrames = floor(sLength/((1-wOverlap)*wLength))-1;

% Compute left and right boundaries
% This is only checked for wOverlap = 1/2, boundary problems may exist
wLeft  = 1:wLength*(1-wOverlap):sLength-wLength+1;
wRight = wLength:wLength*(1-wOverlap):sLength;

% Store all extracted LPC coefficients
allCoeffs = zeros(numFrames, lpcOrder+1);

for i=1:numFrames
    % Pull out current segment and multiply it by window
    curSegment = wav(wLeft(i):wRight(i));
    winSegment = win.*curSegment;

    % The below is different from non-indicator version
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(coastJoint)
        % Compute the energy in the entire signal spectrum.
        energy(i) = sum(winSegment.^2);
    else
        % Compute energy for each energy band 
        psd = abs(fft(winSegment));
        psdRight = psd(1:wLength/2).^2;
        energyLocsInd = floor(energyLocsHz./(fs/2).*length(psdRight));
        for j = 1:size(energyLocsInd,1)
            energy(i,j) = sum(psdRight(energyLocsInd(j,1):energyLocsInd(j,2)));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get the LPC coefficients and associated error
    [allCoeffs(i,:), allError(i)] = lpc(winSegment, lpcOrder);
end

% The below is different from non-indicator version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize indices and compute them
formantInds = ones(numFrames,num_formants);
if coastJoint
    inds = energy < quantile(energy,[quantThresh]);
    formantInds(inds,:) = zeros(sum(inds),num_formants);
else
    for j = 1:numFrames
        inds(j,:) = energy(j,:) < quantile(energy,[quantThresh]);
    end
    % Zero out the indices where coasting should occur
    formantInds(inds) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert LPC Coefficients to cepstral coefficients
[C] = lpc2c(-allCoeffs(:,2:end)',cepOrder);
