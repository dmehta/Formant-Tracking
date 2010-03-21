% Loads the data structure containing both TIMIT and VTR Formant
% information, iterates through all the examples and writes out the
% .wav files to the target directory
%
% INPUT
%   dataSource - Location .mat file with the data
%   targetDir  - directory to which to write out the .wav files
%
% Created on: 03/23/07
% Created by: Daniel Rudoy

function [] = writeTimitDataWav(dataSource, targetDir)

load(dataSource);
fs = 16000; % Sampling rate of TIMIT

for i = 1:length(DATA)
    curData = DATA{i};    
    fileName = strcat(targetDir, '\', 'Timit',int2str(i),'.wav');
   
    % Write out at 16 bits
    wavwrite(curData.timitData, fs, fileName);
end