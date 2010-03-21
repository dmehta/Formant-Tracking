% This code loads the VTR data and corresponding TIMIT segments
% into a single data structure for post processing by formant
% tracker testing software
function [] = generateVtrandTimitExamples()
 
% Location of VTR Formants directory
vtrDataBase = 'C:\Documents and Settings\Dan\My Documents\Research\Speech\Databases\VTRFormants\Train\';
% Location of TIMIT directory
timit_dir = 'C:\Documents and Settings\Dan\My Documents\Research\Speech\Databases\TIMIT\TIMIT\TRAIN';



% Find dialect folders
dialectFolders = dir(vtrDataBase);

kk = 1;
numCurPhone = 1;

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


            file_dir = strcat(timit_dir,'\', curDialect, '\', curSpeaker);
            file_name=strcat(file_dir,'\',curTimitUtter,'.wav');
            
            % Read in sentence
            entireSeg = TIMITread(file_dir,file_name);
            
            % Pull out the content of the sentence
            info_fileName = strcat(file_dir,'\',curTimitUtter,'.txt');
            fid = fopen(info_fileName,'r');
            curLine = fgetl(fid);
            fclose(fid);
            
            % Read in the location of the silence at beginning of sentence
            endInd = find(curTimitUtter == '.');
            curTimitUtter2 = curTimitUtter(1:endInd-1);
            [ind_begSilence,ind_endSilence]=extract_ind_phn(file_dir,curTimitUtter, 'h#');
            
            % Store data into structure for later use by other scripts
            curData.vtrData   = curDataSet;
            curData.timitData = entireSeg; 
            curData.fileName = file_name;
            curData.sentence = curLine;
            curData.initSilence = [ind_begSilence, ind_endSilence];
            DATA{numCurPhone} = curData;
            numCurPhone = numCurPhone + 1; 
        end
    end
end

save VTRandTIMITtrain DATA;