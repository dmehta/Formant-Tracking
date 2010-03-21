% function [state_progression,BW_data] = getDatabaseSeq(filename,num_formants,sample_index,BW_flag)
% Generates Data for testing Deng Outline, for fixed bandwidths, and
% computes the noise to include for a given SNR using data base data
% Also adds the desired noise to the cepstrum coeffs
%
% INPUT
%
% fileName    - url of file that holds all vtr data
% numFormants - number of formants
% sampleIndex - An integer between 1 and 516, indicating which VTR utterance to use
%
% Author: Daniel Spendley, Daniel Rudoy
% Created : 03/3/2007
% Last Updated : 12/13/2007, 03/13/2007


function [f, bw] = vtrFormantRead(filename, numFormants, sampleIndex)

load(filename); % Load all VTR data

vtrData = DATA{sampleIndex}.vtrData'; % Get utterance

% Pull out formants and bandwidths, need to multiply by constants since
% the VTR data is in kHz
f  = 1000*vtrData(1:numFormants, :);  
bw = 500*vtrData(5:(4+numFormants),:); % Scaling on BW somehow a little different