%function [f,bw] = wavesurferFormantRead(fname,num_formants);
% Function to read Wavesurfer 1.8.5 Formant Tracker Output (.frm files)
% 1. Load sound into Wavesurfer; choose "Speech Analysis"
% 2. However over spectrogram and right-click
% 3. Choose "Properties" and adjust formant analyis algorithm inputs
% (no way to tabulate these automatically--even if "header" selected)
% 4. Right-click as in steps 2/3, choose "Save Data File"

%
% Author: Daniel Spendley, Daniel Rudoy
% Created:  03/08/2007
% Modified: 12/13/2007

function [f,bw] = wavesurferFormantRead(fname,numFormants)

% Open file and scan text data for numbers
data = load(fname);
f = data(:,1:numFormants)';
%bw = data(:,(numFormants+1):end)';
bw = data(:,5:1:5+numFormants-1)';
