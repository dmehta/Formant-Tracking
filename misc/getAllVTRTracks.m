% Cycle through the entire VTR database and write contents to a .MAT file
%
% Author:   Daniel Rudoy
% Created:  12/13/2007
% Modified: 12/13/2007

function [] = getAllVTRTracks()

vtrTrain = 'C:\Documents and Settings\Dan\My Documents\Research\Speech\Databases\VTRFormants\Train\';
vtrTest  = 'C:\Documents and Settings\Dan\My Documents\Research\Speech\Databases\VTRFormants\Test\';

sets = {vtrTrain, vtrTest};

numCurPhone = 1;

for d = 1:length(sets)

    % Location of VTR Formants directory (toggle between train and test)
    vtrDataBase = sets{d};

    % Find dialect folders
    dialectFolders = dir(vtrDataBase);

    kk = 1;

    for i = 3:length(dialectFolders)
        curDialect = dialectFolders(i).name;
        curDialectDir = strcat(vtrDataBase,curDialect);
        speakerFolders = dir(curDialectDir);
        for j = 3:length(speakerFolders)

            curSpeaker = speakerFolders(j).name;
            curGender = curSpeaker(1);
            curSpeakerDir = strcat(curDialectDir, '\', curSpeaker);
            utterances = ls(strcat(curSpeakerDir,'\*.fb'));

            [r c] = size(utterances);
            for k = 1:r
                clear vtrFileName;
                vtrFileName = utterances(k,:);
                endInd = find(utterances(k,:) == '.');
                curTimitUtter = vtrFileName(1:endInd-1);
                curUtterance = strcat(curSpeakerDir, '\', vtrFileName);
                curDataSet = loadVTRDatabase(curUtterance);
                dataSets{kk} = curDataSet;
                kk = kk+1;

                % Store data into structure for later use by other scripts
                curData.vtrData   = curDataSet;            
                DATA{numCurPhone} = curData;
                
                numCurPhone = numCurPhone + 1;
            end
        end
    end
end
save ../data/VTR_Timit/allVTRTracks DATA;