%% after running runVTRReplicate.m or loading .mat file

%% A basic plotting routine to visualize results
plotStateTracks(trueStateVTR,estTracks,titleCell);

%% additional plotting routines Nov 2010
figure, plot(trueStateVTR')
hold on, plot(bwDataVTR')
% figure, plot(trueStateWS')
figure, plot(estTracks')

%% state tracks with uncertainties and truth
plotStateTracksFZ_EstVar_Truth(trueStateVTR,estTracks,estVar,titleCell,numFormants,trackBW)
disp(['Mean RMSE: ', num2str(mean(rmse))])
rmse

%% Super-impose Ground Truth over a spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2BW(wav, trueStateVTR, aParams, numAntiF, 0);
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Super-impose EKS over a spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

figure
plotSpecTracks2BW(wav, estTracks, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Wavesurfer tracks on spectrogram
figure
plotSpecTracksWS(wav, trueStateWS, bwDataWS, aParams, trackBW)
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% set(gca, 'PlotBoxAspectRatio', [10 1 1])
 
%% Praat tracks on spectrogram
figure
praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
params = getFormantParamsDefault;
 params.numTracks = numFormants;
 params.winLen = aParams.wLengthMS/1000;
 params.timestep = params.winLen*aParams.wOverlap;
 params.maxformantHz = fs/2;
 params.peHz = params.maxformantHz;
[f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(wav, aParams.fs, praat_dir, params);
plotSpecTracksPraat(wav, f', bw', aParams, formantTimes, trackBW);
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and EKS on spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

figure
plotSpecTracks2BW_Truth(wav, trueStateVTR, estTracks, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and WaveSurfer on spectrogram
figure
plotSpecTracksWS_Truth(wav, trueStateWS, bwDataWS, trueStateVTR, aParams, trackBW)
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% Ground truth and Praat tracks on spectrogram
praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
params = getFormantParamsDefault;
 params.numTracks = numFormants;
 params.winLen = aParams.wLengthMS/1000;
 params.timestep = params.winLen*(1-aParams.wOverlap);
 params.maxformantHz = fs/2;
 params.peHz = params.maxformantHz;
[f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(wav, aParams.fs, praat_dir, params);
plotSpecTracksPraat_Truth(wav, f', bw', trueStateVTR, aParams, formantTimes, trackBW);
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% set(gca, 'PlotBoxAspectRatio', [10 1 1])

%% EKS with uncertainty on spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

figure
plotSpecTracks2_estVar(wav, estTracks, x_errVar, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
% set(gca, 'PlotBoxAspectRatio', [3 1 1])

%% Ground truth and EKS with uncertainty on spectrogram
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

plotSpecTracks2_estVar_Truth(wav, estTracks, x_errVar, trueStateVTR, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))
set(gca, 'PlotBoxAspectRatio', [7 1 1])

%% ...and speech frames
% Sinds(Sinds > length(frameInds)) = []; % added to enforce alignment with frameInds
% xaxe = seconds(frameInds, 1/(aParams.wLengthMS/1000*aParams.wOverlap));
sLength = round((numFrames+1)*wLength/2);
xaxe = 1/fs*(wLength*(1-wOverlap):wLength*(1-wOverlap):sLength-wLength*(1-wOverlap)+1);
xaxeS = xaxe(Sinds);
yaxeS = frameInds(Sinds);
plot(xaxeS, 0, 'g.')

%% ...and phone class color coding
% TIMIT phone classes: 
% silence (0, none), vowel (1, blue), semivowel (2, green), nasal (3, cyan), fricative (4, magenta),
% affricate (5, red), stop (6, black)
sLength = round((numFrames+1)*wLength/2);
xaxe = 1/fs*(wLength*(1-wOverlap):wLength*(1-wOverlap):sLength-wLength*(1-wOverlap)+1);
classColor = {'b', 'g', 'c', 'm', 'r', 'k'};
for class = 1:length(classColor)
    Cinds = frameInds == class;
    xaxeC = xaxe(Cinds);
    if ~isempty(xaxeC) % only plot if class exists
        plot(xaxeC, 0, [classColor{class}, '.'])
    end
end

%% Ground truth and EKS on spectrogram...and speech frames
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

figure
plotSpecTracks2BW_Truth(wav, trueStateVTR, estTracks, aParams, numAntiF, trackBW);
% plotSpecTracks2BW(wav, estTracks, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))

Sinds(Sinds > length(frameInds)) = []; % added to enforce alignment with frameInds
sLength = round((numFrames+1)*wLength/2);
xaxe = 1/fs*(wLength*(1-wOverlap):wLength*(1-wOverlap):sLength-wLength*(1-wOverlap)+1);
xaxeS = xaxe(Sinds);
yaxeS = frameInds(Sinds);
plot(xaxeS, 0, 'g.')
% ylim([0 5000])

%% Ground truth and EKS on spectrogram...and phone class color coding
aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

figure
plotSpecTracks2BW_Truth(wav, trueStateVTR, estTracks, aParams, numAntiF, trackBW);
% plotSpecTracks2BW(wav, estTracks, aParams, numAntiF, trackBW);
firstline = get(get(gca,'Title'),'String');
title(strrep([{['File: ', dataFileNameIn]}; firstline], '_', '\_'))

Sinds(Sinds > length(frameInds)) = []; % added to enforce alignment with frameInds
sLength = round((numFrames+1)*wLength/2);
xaxe = 1/fs*(wLength*(1-wOverlap):wLength*(1-wOverlap):sLength-wLength*(1-wOverlap)+1);
xaxeS = xaxe(Sinds);
yaxeS = frameInds(Sinds);
classColor = {'b', 'g', 'c', 'm', 'r', 'k'};
for class = 1:length(classColor)
    Cinds = frameInds == class;
    xaxeC = xaxe(Cinds);
    if ~isempty(xaxeC) % only plot if class exists
        plot(xaxeC, 0, [classColor{class}, '.'])
    end
end