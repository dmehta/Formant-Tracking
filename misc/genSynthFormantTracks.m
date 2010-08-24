function [frmVals, bwVals] = genSynthFormantTracks(pNoiseVar, numObser, frmVals0, bwVals0)

% Starting from initial formant states generate a synthetic formant track
% via a Normal random walk with variance given by pNoiseVar; N.B. variance of bandwidth
% tracks scaled by 1/25

% Author: Daniel Rudoy 
% Created : March 2007
% Modified: 12/13/2010


numStates = length(frmVals0); % Number of formants tracks to generate

frmVals = zeros(numStates, numObser); % Initialize matrix for storing formant tracks
bwVals  = zeros(numStates, numObser); % Initialize matrix for storing formant tracks

frmVals(:,1) = frmVals0'; % Initial formant values
bwVals(:,1)  = bwVals0';  % Initial bandwidth values

% Evolve formants and bandwidths (keep variance of bandwidths very low)
for i = 2:numObser
    frmVals(:,i) = abs(frmVals(:,i-1)) + sqrt(pNoiseVar)*randn(numStates,1);
    bwVals(:,i) = abs(bwVals(:,i-1)) + sqrt(pNoiseVar/25)*randn(numStates,1);
end