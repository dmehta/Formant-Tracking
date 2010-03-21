% Test script to see if we can load data from VTR database

% This code loads only the VTR Labels does not load the original signal
% The corrected data includes only F1, F2, F3 -- not F4, and B1-B4

function [dataSet] = loadVTRDatabase(fileName)

fid=fopen(fileName, 'r', 'b');
n_frame=fread(fid, 1, 'int32');
samPeriod=fread(fid, 1, 'int32');
sampSize=fread(fid, 1, 'int16');
numComps=sampSize/4;
fileType=fread(fid, 1, 'int16');
dataSet=zeros(n_frame,numComps);
for n = 1:n_frame
    a=fread(fid, numComps, 'float');
    dataSet(n,:)=a';
end
fclose(fid);