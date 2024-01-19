% Phonetic classification of each frame using TIMIT phone labels. Frame is
% classified as specific phone class if majority of samples within frame is
% labeled as a given phone. If multiple phones tie, priority will be given
% in the following order: vowel, semivowel, nasal, fricative, affricate,
% stop.
% 
% Silent phone labels per sample:
%   Silence: pau (pause), epi, (epenthetic pause), h# (silence at beginning and end of utterance)
%   Stop closure: bcl, dcl, gcl, pcl, tcl, kcl
%   Glottal stop: q
%
% INPUT
% wav            - waveform
% cellPhones     - cell array containing all phonemes in utterance
% analysisParams - Datastructure specifying the short-time analysis scheme
%                  within which to do detection
% doPlot         - Show the mask that is produced
%
% OUTPUT
% frameInds      - vector the length of the number of frames, with value at each frame 
%                   indicating one of the following phonetic
%                   classifications: vowel (1), semivowel (2), nasal (3),
%                   fricative (4), affricate (5), stop (6), silence (0)
% 
% Author: Daryush Mehta
% Created:  11/21/2010
% Modified: 10/17/2011

function frameInds = TIMITphnclassify(wav, cellPhones, analysisParams, doPlot)
fs_in = 16e3;

% Extract the necessary parameters from the analysis set
wLength           = analysisParams.wLength;    % Window length
wOverlap          = analysisParams.wOverlap;   % Window overlap
fs                = analysisParams.fs;         % Sampling rate

% Truncate samples that don't fall within a fixed number of windows
% (PJW: Seems this is what VTR database/Wavesurfer does)
%wav = wav(1:end-mod(length(wav),wLength)); %DGR: this might be an issue
numFrames = floor(length(wav)/(wLength/2))-1;
sLength = round((numFrames+1)*wLength/2);

% Compute number of frames
% numFrames = floor(sLength/((1-wOverlap)*wLength))-1;

% Compute left and right boundaries
% This is only checked for wOverlap = 1/2, boundary problems may exist
wLeft  = 1:wLength*(1-wOverlap):sLength-wLength+1;
wRight = wLength:wLength*(1-wOverlap):sLength;

% Pre-allocate array
frameInds = zeros(1,numFrames);

% 6 phonetic classes and 1 silence class
vowelphn = {'iy', 'ih', 'eh', 'ey', 'ae', 'aa', ...
    'aw', 'ay', 'ah', 'ao', 'oy', 'ow', 'uh', ...
    'uw', 'ux', 'er', 'ax', 'ix', 'axr'};
semivowelphn = {'l', 'r', ...
    'w', 'y', ...
    'hh', 'hv', ...
    'el'};
nasalphn = {'m', 'n', 'ng', ...
    'em', 'en', 'eng', ...
    'nx'};
fricphn = {'s', 'sh', 'z', 'zh', 'f', 'th', 'v', 'dh'};
affricphn = {'jh', 'ch'};
stopphn = {'b', 'd', 'g', 'p', 't', 'k', 'dx'};
silencephn = {'pau', 'epi', 'h#', ...
    'bcl', 'dcl', 'gcl', 'pcl', 'tcl', 'kcl', ...
    'q'};

for i=1:numFrames
    % Count number of silence samples within this frame
    count = [0 0 0 0 0 0 0]; % order: vowel, semivowel, nasal, fricative, affricate, stop, silence
    for jj = wLeft(i):wRight(i)
        lookupInd = floor(jj*fs_in/fs);
        if lookupInd <= length(cellPhones) % NB phone labeling might be short per timit readme, set rest to 0      
            if any(strcmpi(cellPhones{lookupInd}, vowelphn))
                count(1) = count(1) + 1;
            end
            if any(strcmpi(cellPhones{lookupInd}, semivowelphn))
                count(2) = count(2) + 1;
            end
            if any(strcmpi(cellPhones{lookupInd}, nasalphn))
                count(3) = count(3) + 1;
            end
            if any(strcmpi(cellPhones{lookupInd}, fricphn))
                count(4) = count(4) + 1;
            end
            if any(strcmpi(cellPhones{lookupInd}, affricphn))
                count(5) = count(5) + 1;
            end
            if any(strcmpi(cellPhones{lookupInd}, stopphn))
                count(6) = count(6) + 1;
            end
            if any(strcmpi(cellPhones{lookupInd}, silencephn))
                count(7) = count(7) + 1;
            end
        else
            count(7) = count(7) + 1; % this sets samples outside timit phone labeling to silence
        end 
    end
    if count(7) < wLength - 1 % then the frame is a speech frame, else frame value is 0
        [val, index] = max(count(1:end-1));
        frameInds(i) = index(1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(doPlot)
    figure;
    imagesc(frameInds);
    title('Result of TIMIT Speech Classification');
    xlabel('Short-time frame number');
end