% getDefault.m
%
% Author: Daryush Mehta
% Created: 07/29/2011
% Last modified: 08/01/2011
% 
% Helper function to load default parameters for KARMA approach to formant
% and antiformant tracking
% 
% function [numFormants, numAntiF, aParams, cepOrder, cepType, algFlag, x0, bwStates, formantInds, trackBW] = getDefault()
% 
% INPUT:
%    (none)
% 
% OUTPUT:
%    numFormants:   Number of formants to track; default: 3
%    numAntiF:      Number of anti-formants to track; default: 0
%    aParams:       Analysis parameters, structure with one or more of the following fields:
%                           peCoeff     % Pre-emphasis coefficient; default: 0.7
%                           wType       % window type; default: 'Hamming'
%                           wLengthMS   % Length of window, in ms; default: 20
%                           wOverlap    % Window overlap factor; default: 0.5
%                           lpcOrder    % Number of AR coefficients; default: 12
%                           zOrder      % Number of MA coefficients; default: 0
%                           fs          % Downsampling frequency, in Hz; default: 7000
%    cepOrder:      Number of cepstral coefficients to observe; default: 15
%    cepType:       Type of cepstrum to compute; 1 = ARMA, 2 = real
%                           default: 1
%    algFlag:       Flag for algorithm to run, 1 for EKS, 2 for EKF;
%                           default: 1
%    x0:            Column vector of initialized center frequency values
%                           for formants then antiformants, in Hz
%                           Length is (numFormants + numAntiformants);
%                           default: [500 1500 2500 ... ] for formants
%                           [1000 2000 3000 ... ] for antiformants
%    bwStates:      Column vector of initialized bandwidth values
%                           for formants then antiformants, in Hz
%                           Length is (numFormants + numAntiformants);
%                           default: [80 120 160 ... ] for formants
%                           [80 120 160 ... ] for antiformants
%    formantInds:   Use speech activity detection? 1 for yes (tracks will be
%                           coasted during non-speech frames), [] for no; 
%                           can also input masking matrix of
%                           formant indices: 0 if masked, 1 if
%                           unmasked; default: 1
%    trackBW:       Flag to calculate bandwidths (1) or fix to
%                           initial values in bwStates (2)
%                           default: 1

function [numFormants, numAntiF, aParams, cepOrder, cepType, algFlag, x0, bwStates, formantInds, trackBW] = getDefault()

numFormants = 3;
numAntiF = 0;
aParams.peCoeff = 0.7;
aParams.wType = 'hamming';
aParams.wLengthMS = 20;
aParams.wOverlap = 0.5;
aParams.lpcOrder = 12;
aParams.zOrder = 0;
aParams.fs = 7000;
cepOrder = 15;
cepType = 1;
algFlag = 2;

x0 = 500 + 1000*(0:(numFormants - 1))';
x02 = 1000 + 1000*(0:(numAntiF - 1))';
x0 = [x0; x02];

bwStates = 80 + 40*(0:(numFormants - 1))';
bwStates2 = 80 + 40*(0:(numAntiF - 1))';
bwStates = [bwStates; bwStates2];

formantInds = 1;
trackBW = 1;