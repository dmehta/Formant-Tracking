function run_make_VTR_Results_All()
%function genAllVTRFigs()
% Generate .pdf files for all utterances in VTR database,
% according to settings used in (Mehta, Rudoy, Wolfe, 2012).
% 
% Plots are labeled with utterance and other TIMIT info (gender, dialog,
% etc), and TIMIT phoneme segmentation boundaries are indicated.
%
% Make PDF displaying all 516 VTR results

% if ~exist('plotSpecTracks2BW_Truth', 'file')
%     addpath(genpath('../')); % executes if directories not on path
% end

% % VTR, source 1:
% root_fname = 'VTR';
% savedir = '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';
% savedir = '..\results\EKS_WS_Praat_trackBW0_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';

% VTRsynth3500, source 4:
% root_fname = 'VTRsynth3500';
% savedir = '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';

% VTRsynthf03500, source 5:
% root_fname = 'VTRsynthf03500';
% savedir = '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';

% set up loop names
sources = [1, 1, 4, 5];
root_fnames = {'VTR', 'VTR', 'VTRsynth3500', 'VTRsynthf03500'};
savedirs = { ...
    '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix', ...
    '..\results\EKS_WS_Praat_trackBW0_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix', ...
    '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix', ...
    '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix'...
    };

for loop = 1:length(root_fnames)
    root_fname = root_fnames{loop};
    savedir = savedirs{loop};
    
    for ii = 1:516
        figure(1)
        load(fullfile(savedir, [root_fname, num2str(ii)]))

        %% 
        subplot(411)
        plotSpecBlank(wav, aParams);
        % display TIMIT sentence information
        covs = getTimitCovariates(ii);
        title(strrep({dataFileNameIn; [covs.dialect ' ' covs.gender ': "' covs.sentence '"']}, '_', '\_'))

        %% Ground truth and EKS on spectrogram...and speech frames
        subplot(412)
        aParams.wLength = floor(aParams.wLengthMS/1000*aParams.fs);
        aParams.wLength = aParams.wLength + (mod(aParams.wLength,2)~=0); % Force even

        plotSpecTracks2BW_Truth(wav, trueStateVTR, estTracks, aParams, numAntiF, trackBW);
        firstline = get(get(gca,'Title'),'String');
        title([{['EKS: RMSE = ', num2str(round(rmseAllS(ii,1,1))), ', ', ...
            num2str(round(rmseAllS(ii,1,2))), ', ', ...
            num2str(round(rmseAllS(ii,1,3))), ', ', ...
            num2str(round(rmseAllS(ii,1,4))), ...
            ' Hz', ...
            ]}; firstline])

        %% Ground truth and WaveSurfer on spectrogram
        subplot(413)
        plotSpecTracksWS_Truth(wav, trueStateWS, bwDataWS, trueStateVTR, aParams, trackBW)
        firstline = get(get(gca,'Title'),'String');
        title([{['WaveSurfer: RMSE = ', num2str(round(rmseAllS(ii,2,1))), ', ', ...
            num2str(round(rmseAllS(ii,2,2))), ', ', ...
            num2str(round(rmseAllS(ii,2,3))), ', ', ...
            num2str(round(rmseAllS(ii,2,4))), ...
            ' Hz', ...
            ]}; firstline])

        %% Ground truth and Praat tracks on spectrogram
        subplot(414)
        praat_dir = 'C:\Documents and Settings\Daryush\My Documents\MATLAB';
        %praat_dir = 'Y:\matlab\chsv';
        params = getFormantParamsDefault;
         params.numTracks = numFormants;
         params.winLen = aParams.wLengthMS/1000;
         params.timestep = params.winLen*(1-aParams.wOverlap);
         params.maxformantHz = fs/2;
         params.peHz = params.maxformantHz;
         if (sources(loop) == 1 && (vtrDbNum == 153 || vtrDbNum == 193)) || ...
             (sources(loop) == 3 && (vtrDbNum == 78 || vtrDbNum == 464)) || ...
             (sources(loop) == 5 && (vtrDbNum == 124 || vtrDbNum == 131 || vtrDbNum == 138 || vtrDbNum == 217 || vtrDbNum == 228 || vtrDbNum == 416))
             params.method = 'sl';
         end
        [f, bw, numberOfFrames, frameSpacingTime, startFrameTime, formantTimes] = getFormantTrack(wav, aParams.fs, praat_dir, params);
        plotSpecTracksPraat_Truth(wav, f', bw', trueStateVTR, aParams, formantTimes, trackBW);
        firstline = get(get(gca,'Title'),'String');
            title([{['Praat: RMSE = ', num2str(round(rmseAllS(ii,3,1))), ', ', ...
            num2str(round(rmseAllS(ii,3,2))), ', ', ...
            num2str(round(rmseAllS(ii,3,3))), ', ', ...
            num2str(round(rmseAllS(ii,3,4))), ...
            ' Hz', ...
            ]}; firstline])

        %% Add stuff to each plot
        for plotnum = 1:4
            subplot(4,1,plotnum)
            % aspect ratio
            set(gca, 'PlotBoxAspectRatio', [10 1 1])

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
        end

        subplot(4,1,4)
        text(0, -5000, {'RMSE is over speech-labeled frames for F1, F2, F3, and Average, relative to VTR database truth. '; ...
            'Frames color coded by TIMIT phone class: vowel (blue), semivowel (green), nasal (cyan), '; ...
            'fricative (magenta), affricate (red), stop (black)'})

        %% print to PDF
        set(gcf,'PaperPosition',[0 0 8.5 11]);

        %set(gcf, 'Units', 'Normalized', 'Position', [0 0 1 1])
        %print(1, '-djpeg', fullfile(savedir, [root_fname, num2str(ii)]))
        % Print to pdf using default resolution (684--'-r864')
        print('-dpdf', fullfile(savedir, [root_fname, num2str(ii), '.pdf']));
        
        % now print to postscript so we can append pages; doesn't handle
        % large files well
        %if ii == 1
        %    print ('-dpsc2', fullfile(savedir, [root_fname, '.ps']))
        %else
        %    print ('-dpsc2', fullfile(savedir, [root_fname, '.ps']), '-append')
        %end

        close(1)
    end
    
    % convert postscript file to PDF
    % http://www.mathworks.com/matlabcentral/fileexchange/19516-ps2pdf
    %     ps2pdf('psfile', fullfile(savedir, [root_fname, '.ps']), 'pdffile', ...
    %         fullfile(savedir, [root_fname, '_NEW.pdf']), ...
    %         'gscommand', 'C:\Program Files\gs\gs9.02\bin\gswin32c.exe', ...
    %         'gslibpath', 'C:\Program Files\gs\gs9.02\lib', ...
    %         'gsfontpath', 'C:\Program Files\gs\gs9.02\lib', ...
    %         'gspapersize', 'letter', ...
    %         'deletepsfile', 1);
end

%%
function timitCovs = getTimitCovariates(timitNum)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load Timit utterance matfile
matFile = ['Timit' int2str(timitNum) '.mat'];
load(matFile);

% Get utterance text
sentenceLetterInd = find(isletter(data.sentence));
sentence = data.sentence(sentenceLetterInd(1):end);

% Get speaker gender
fnameFilesepInd = strfind(data.fileName,filesep);
genderInitial = data.fileName(fnameFilesepInd(end-1)+1);
switch lower(genderInitial)
    case 'f'
        gender = 'female';
    case 'm'
        gender = 'male';
    otherwise
        error('Error parsing gender string: neither male nor female');
end

% Get speaker dialect
dialectCode = data.fileName(fnameFilesepInd(end-1)-1);
switch lower(dialectCode)
    case '1'
        dialect = 'New England';
    case '2'
        dialect = 'Northern';
    case '3'
        dialect = 'North Midland';
    case '4'
        dialect = 'South Midland';
    case '5'
        dialect = 'Southern';
    case '6'
        dialect = 'New York City';
    case '7'
        dialect = 'Western';
    case '8'
        dialect = 'Army Brat';
    otherwise
        error('Error parsing dialect string: not between 1 and 8');
end

% Assign outputs and return
timitCovs.sentence = sentence;
timitCovs.gender = gender;
timitCovs.dialect = dialect;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%