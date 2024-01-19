% Speech activity detection using TIMIT phone labels. Frame is considered
% speech if majority of samples within frame is labeled as a non-silent
% phone.
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
% Author: Daryush Mehta
% Created:  11/21/2010
% Modified: 11/21/2010
function frameInds = timitVAD(wav, cellPhones, analysisParams, doPlot)
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
silencephn = {'pau', 'epi', 'h#', ...
    'bcl', 'dcl', 'gcl', 'pcl', 'tcl', 'kcl', ...
    'q'};

for i=1:numFrames
    % Count number of silence samples within this frame
    count = 0;
    for jj = wLeft(i):wRight(i)
        lookupInd = floor(jj*fs_in/fs);
        if lookupInd <= length(cellPhones) % NB phone labeling might be short per timit readme, set rest to 0
            if any(strcmpi(cellPhones{lookupInd}, silencephn))
                count = count + 1;
            end
        else
            count = count + 1; % set rest of samples to silence, increment count!
        end 
    end
    if count < wLength - 1 % label speech
        frameInds(i) = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(doPlot)
    figure;
    imagesc(frameInds);
    title('Result of TIMIT Speech Activity Detector');
    xlabel('Short-time frame number');
end