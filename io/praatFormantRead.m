%function [f,bw] = praatFormantRead(fname,num_formants);
% Function to read Praat 4.22 Formant Tracker Output (.Formant files)
% 1. Load sound into Praat; choose "Formants & LPC"
% 2. Navigate to "To Formant (hack)"; then "To Formant (sl)"
% Choose default and compute formants
% (Navigate to "Draw" to display formants)
% Navigate to "Write" on main menu; 
% Select "Write to short text file"

% Author: Patrick J. Wolfe, Daniel Rudoy
% Created:  03/08/2007
% Modified: 12/23/2007

function [f, bw] = praatFormantRead(fname, numFormants)

% Open file and scan text data for numbers
fid = fopen(fname);
fdata = textscan(fid,'%f','headerLines',11);
fclose(fid);

% Formants
f = fdata{1}(1:2:end);
f((numFormants+1):(numFormants+1):end) = [];

numFrames = length(f)/numFormants;
if (mod(numFrames,1) ~=0),
    error('Incorrect indexing: non-integer number of frames')
end
f = reshape(f,numFormants,numFrames);

% Bandwidths
bw = fdata{1}(2:2:end);
bw((numFormants+1):(numFormants+1):end) = [];
bw = reshape(bw,numFormants,numFrames);
