%% use after loading VTR516.mat to find error metrics for all, males, and
%% females OR for individual speaker errors, load VTR3.mat (e.g.)

%% Calculate error, used for JASA submission outputs (see RMSEAllS for individual
%% RMSE values)
% % load('..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWws\VTR516.mat')
% indices = 1:516; % for VTRSynth
% % indices = [1:77 79:463 465:516]; % for praat
% % indices = [1:123 125:130 132:216 218:415 417:516]; % for praat on VTRsynthf03500
% 
% % rmseAllS(vtrDbNum, trackerType, numFormants_error+1)
% % 1. EKS, 2. WaveSurfer, 3. Praat
% % indices = 1:516;
% rmse_each_table = zeros(numTrackers, numFormants_error+1);
% for ii = trackersToRun % # trackers
%     for jj = 1:4 % # rmse values
%         tmp = rmseAllS(indices, ii, jj);
%         rmse_each_table(ii, jj) = mean(tmp(~isnan(tmp))); % take into account any NaNs
%     end
% end
% 
% % males are blue, females are red
% figure, stem(indices, rmseAllS(indices, 1, 1), 'o') % F1 rmse
% hold on, stem(indices(logical(Fvect)), rmseAllS(indices(logical(Fvect)), 1, 1), 'ro') % FEMALES
% xlabel('VTR number')
% ylabel('RMSE (Hz)')
% title('RMSE of F1 (Hz), Male (blue) and female (red)')
% 
% figure, stem(indices, rmseAllS(indices, 1, 2), 'o') % F2 rmse
% hold on, stem(indices(logical(Fvect)), rmseAllS(indices(logical(Fvect)), 1, 2), 'ro') % FEMALES
% xlabel('VTR number')
% ylabel('RMSE (Hz)')
% title('RMSE of F2 (Hz), Male (blue) and female (red)')
% 
% figure, stem(indices, rmseAllS(indices, 1, 3), 'o') % F3 rmse
% hold on, stem(indices(logical(Fvect)), rmseAllS(indices(logical(Fvect)), 1, 3), 'ro') % FEMALES
% xlabel('VTR number')
% ylabel('RMSE (Hz)')
% title('RMSE of F3 (Hz), Male (blue) and female (red)')
% 
% figure, stem(indices, rmseAllS(indices, 1, 4), 'o') % overall rmse
% hold on, stem(indices(logical(Fvect)), rmseAllS(indices(logical(Fvect)), 1, 4), 'ro') % FEMALES
% xlabel('VTR number')
% ylabel('RMSE (Hz)')
% title('Average RMSE over 3 formants, Male (blue) and female (red)')
% 
% rmse_each_table'
% round(rmse_each_table')

%% cycle through VTR database and save gender (female flag) to
%% TIMITgender.mat (only need to do once before executing code below)
% vtrDbNums = 1:516;
% Fvect = zeros(1, length(vtrDbNums));
% 
% for vtrDbNum = vtrDbNums
%     matFileName = strcat('../data/VTR_Timit/Timit',num2str(vtrDbNum),'.mat');
%     load(matFileName)
%     [x, y, z] = fileparts(data.fileName);
%     gender = x(end-4);
%     if strcmpi(gender, 'f')
%         Fvect(vtrDbNum) = 1;
%     end
% end
% 
% save(strcat('../data/VTR_Timit/TIMITgender.mat'), 'Fvect')

        %% 1. For JASA re-submission: Calculate error within phonetic classes
        % loading _error.mat file replaces this cell, or load VTR516.mat and re-run

        trueStateVTR2_All_Db = cell(numTrackers, numFormants_error, numClasses);
        estTracks2_All_Db = cell(numTrackers, numFormants_error, numClasses);
        trueStateVTR2_All_Db_class = cell(numTrackers, numFormants_error);
        estTracks2_All_Db_class = cell(numTrackers, numFormants_error);
        trueStateVTR2_All_Db_class_f = cell(1, numTrackers);
        estTracks2_All_Db_class_f = cell(1, numTrackers);
        for trackerType = trackersToRun
            for j = 1:numFormants_error
                for class = 1:numClasses
                    for vtrDbNum = vtrDbNums % *** set to vtrDbNum if only interested in current utterance RMSE (otherwise, vtrDbNums)
                    %for vtrDbNum = [1:329 331:516] % VTR F3 blatantly incorrect so skews results
                        % reshape so that we collapse tracks across all utterances, preserving
                        % categories of formant number and phoneme class (and
                        % tracker type)
                        trueStateVTR2_All_Db{trackerType,j,class} = [trueStateVTR2_All_Db{trackerType,j,class} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                        estTracks2_All_Db{trackerType,j,class}    = [estTracks2_All_Db{trackerType,j,class}    estTracks2_All{vtrDbNum,trackerType,j,class}];

                        % reshape so that we collapse tracks across all utterances and all phonetic classes, preserving
                        % category of formant number (and tracker type)                
                        trueStateVTR2_All_Db_class{trackerType,j} = [trueStateVTR2_All_Db_class{trackerType,j} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                        estTracks2_All_Db_class{trackerType,j}    = [estTracks2_All_Db_class{trackerType,j}    estTracks2_All{vtrDbNum,trackerType,j,class}];                

                        % reshape so that we collapse tracks across all utterances, all phonetic classes, and all formants, preserving
                        % tracker type
                        trueStateVTR2_All_Db_class_f{trackerType} = [trueStateVTR2_All_Db_class_f{trackerType} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                        estTracks2_All_Db_class_f{trackerType}    = [estTracks2_All_Db_class_f{trackerType}    estTracks2_All{vtrDbNum,trackerType,j,class}];                
                    end
                end
            end
        end

        % Run this cell after loading VTR_error.mat, VTR_errormale.mat, VTR_errorfemale.mat or running one of previous cells 
        % Number of frames for each phonetic class
        for ii = 1:6, temp = trueStateVTR2_All_Db{1,1,ii};disp(length(temp)), end

        % Number of frames for each formant (should be the same)
        for ii = 1:3, temp = trueStateVTR2_All_Db_class{1,ii};disp(length(temp)), end

        % Number of frames for each tracker
        for ii = 1:3, temp = trueStateVTR2_All_Db_class_f{ii};disp(length(temp)), end

        % now, find error for each formant trajectory and put in table a la Deng (2006)
        % Row 1 (Class 1): F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
        % Row 2 (Class 2): F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
        % ...
        rmse_table = zeros(numClasses, numTrackers*numFormants_error);
        for class = 1:numClasses
            for trackerType = trackersToRun
                for j = 1:numFormants_error
                    differ = estTracks2_All_Db{trackerType,j,class} - trueStateVTR2_All_Db{trackerType,j,class};
                    differ = differ(~isnan(differ)); % get rid of NaNs
                    % calculate RMSE
                    rmse_table(class, j+(trackerType-1)*numTrackers) = sqrt(mean((differ).^2));

                    % calculate mean absolute error
                    %rmse_table(class, j+(trackerType-1)*numTrackers) = mean(abs(differ));
                end
            end
        end
        round(rmse_table)

        % now, find error for each formant trajectory and put in table
        % F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
        rmse_table2 = zeros(1, numTrackers*numFormants_error);
        for trackerType = trackersToRun
            for j = 1:numFormants_error
                differ = estTracks2_All_Db_class{trackerType,j} - trueStateVTR2_All_Db_class{trackerType,j};
                differ = differ(~isnan(differ)); % get rid of NaNs
                % calculate RMSE
                rmse_table2(j+(trackerType-1)*numTrackers) = sqrt(mean((differ).^2));

                % calculate mean absolute error
                %rmse_table2(j+(trackerType-1)*numTrackers) = mean(abs(differ));
            end
        end
        round(rmse_table2)

        % now, find overall error for each tracker and put in table
        % F1, F2, F3 (for all trackers)
        rmse_table3 = zeros(1, numTrackers);
        for trackerType = trackersToRun
            differ = estTracks2_All_Db_class_f{trackerType} - trueStateVTR2_All_Db_class_f{trackerType};
            differ = differ(~isnan(differ)); % get rid of NaNs
            % calculate RMSE
            rmse_table3(trackerType) = sqrt(mean((differ).^2));

            % calculate mean absolute error
            %rmse_table3(trackerType) = mean(abs(differ));
        end
        round(rmse_table3)

        % pick appropriate save function
        save(fullfile(savedir, [prefix, '_error']))
        % save(fullfile(savedir, [prefix, '_errormale']))
        % save(fullfile(savedir, [prefix, '_errorfemale']))


        %% 2. For JASA re-submission: Calculate error within phonetic classes MALES
        % loading _errormale.mat file replaces this cell

        load(strcat('../data/VTR_Timit/TIMITgender.mat'), 'Fvect')

        trueStateVTR2_All_Db = cell(numTrackers, numFormants_error, numClasses);
        estTracks2_All_Db = cell(numTrackers, numFormants_error, numClasses);
        trueStateVTR2_All_Db_class = cell(numTrackers, numFormants_error);
        estTracks2_All_Db_class = cell(numTrackers, numFormants_error);
        trueStateVTR2_All_Db_class_f = cell(1, numTrackers);
        estTracks2_All_Db_class_f = cell(1, numTrackers);
        for trackerType = trackersToRun
            for j = 1:numFormants_error
                for class = 1:numClasses
                    for vtrDbNum = vtrDbNums % *** set to vtrDbNum if only interested in current utterance RMSE (otherwise, vtrDbNums)
                    %for vtrDbNum = [1:329 331:516] % VTR F3 blatantly incorrect so skews results    
                        % only add to master vectors if correct gender
                        if ~Fvect(vtrDbNum) % get only male values
                            % reshape so that we collapse tracks across all utterances, preserving
                            % categories of formant number and phoneme class (and
                            % tracker type)
                            trueStateVTR2_All_Db{trackerType,j,class} = [trueStateVTR2_All_Db{trackerType,j,class} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                            estTracks2_All_Db{trackerType,j,class}    = [estTracks2_All_Db{trackerType,j,class}    estTracks2_All{vtrDbNum,trackerType,j,class}];

                            % reshape so that we collapse tracks across all utterances and all phonetic classes, preserving
                            % category of formant number (and tracker type)                
                            trueStateVTR2_All_Db_class{trackerType,j} = [trueStateVTR2_All_Db_class{trackerType,j} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                            estTracks2_All_Db_class{trackerType,j}    = [estTracks2_All_Db_class{trackerType,j}    estTracks2_All{vtrDbNum,trackerType,j,class}];                

                            % reshape so that we collapse tracks across all utterances, all phonetic classes, and all formants, preserving
                            % tracker type
                            trueStateVTR2_All_Db_class_f{trackerType} = [trueStateVTR2_All_Db_class_f{trackerType} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                            estTracks2_All_Db_class_f{trackerType}    = [estTracks2_All_Db_class_f{trackerType}    estTracks2_All{vtrDbNum,trackerType,j,class}];
                        end
                    end
                end
            end
        end

        %% Run this cell after loading VTR_error.mat, VTR_errormale.mat, VTR_errorfemale.mat or running one of previous cells 
        % Number of frames for each phonetic class
        for ii = 1:6, temp = trueStateVTR2_All_Db{1,1,ii};disp(length(temp)), end

        % Number of frames for each formant (should be the same)
        for ii = 1:3, temp = trueStateVTR2_All_Db_class{1,ii};disp(length(temp)), end

        % Number of frames for each tracker
        for ii = 1:3, temp = trueStateVTR2_All_Db_class_f{ii};disp(length(temp)), end

        % now, find error for each formant trajectory and put in table a la Deng (2006)
        % Row 1 (Class 1): F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
        % Row 2 (Class 2): F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
        % ...
        rmse_table = zeros(numClasses, numTrackers*numFormants_error);
        for class = 1:numClasses
            for trackerType = trackersToRun
                for j = 1:numFormants_error
                    differ = estTracks2_All_Db{trackerType,j,class} - trueStateVTR2_All_Db{trackerType,j,class};
                    differ = differ(~isnan(differ)); % get rid of NaNs
                    % calculate RMSE
                    rmse_table(class, j+(trackerType-1)*numTrackers) = sqrt(mean((differ).^2));

                    % calculate mean absolute error
                    %rmse_table(class, j+(trackerType-1)*numTrackers) = mean(abs(differ));
                end
            end
        end
        round(rmse_table)

        % now, find error for each formant trajectory and put in table
        % F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
        rmse_table2 = zeros(1, numTrackers*numFormants_error);
        for trackerType = trackersToRun
            for j = 1:numFormants_error
                differ = estTracks2_All_Db_class{trackerType,j} - trueStateVTR2_All_Db_class{trackerType,j};
                differ = differ(~isnan(differ)); % get rid of NaNs
                % calculate RMSE
                rmse_table2(j+(trackerType-1)*numTrackers) = sqrt(mean((differ).^2));

                % calculate mean absolute error
                %rmse_table2(j+(trackerType-1)*numTrackers) = mean(abs(differ));
            end
        end
        round(rmse_table2)

        % now, find overall error for each tracker and put in table
        % F1, F2, F3 (for all trackers)
        rmse_table3 = zeros(1, numTrackers);
        for trackerType = trackersToRun
            differ = estTracks2_All_Db_class_f{trackerType} - trueStateVTR2_All_Db_class_f{trackerType};
            differ = differ(~isnan(differ)); % get rid of NaNs
            % calculate RMSE
            rmse_table3(trackerType) = sqrt(mean((differ).^2));

            % calculate mean absolute error
            %rmse_table3(trackerType) = mean(abs(differ));
        end
        round(rmse_table3)

        % pick appropriate save function
        %save(fullfile(savedir, [prefix, '_error']))
        save(fullfile(savedir, [prefix, '_errormale']))
        % save(fullfile(savedir, [prefix, '_errorfemale']))


        %% 3. For JASA re-submission: Calculate error within phonetic classes FEMALES
        % loading _errorfemale.mat file replaces this cell

        load(strcat('../data/VTR_Timit/TIMITgender.mat'), 'Fvect')

        trueStateVTR2_All_Db = cell(numTrackers, numFormants_error, numClasses);
        estTracks2_All_Db = cell(numTrackers, numFormants_error, numClasses);
        trueStateVTR2_All_Db_class = cell(numTrackers, numFormants_error);
        estTracks2_All_Db_class = cell(numTrackers, numFormants_error);
        trueStateVTR2_All_Db_class_f = cell(1, numTrackers);
        estTracks2_All_Db_class_f = cell(1, numTrackers);
        for trackerType = trackersToRun
            for j = 1:numFormants_error
                for class = 1:numClasses
                    for vtrDbNum = vtrDbNums % *** set to vtrDbNum if only interested in current utterance RMSE (otherwise, vtrDbNums)
                    %for vtrDbNum = [1:329 331:516] % VTR F3 blatantly incorrect so skews results    
                        % only add to master vectors if correct gender
                        if Fvect(vtrDbNum) % get only female values
                            % reshape so that we collapse tracks across all utterances, preserving
                            % categories of formant number and phoneme class (and
                            % tracker type)
                            trueStateVTR2_All_Db{trackerType,j,class} = [trueStateVTR2_All_Db{trackerType,j,class} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                            estTracks2_All_Db{trackerType,j,class}    = [estTracks2_All_Db{trackerType,j,class}    estTracks2_All{vtrDbNum,trackerType,j,class}];

                            % reshape so that we collapse tracks across all utterances and all phonetic classes, preserving
                            % category of formant number (and tracker type)                
                            trueStateVTR2_All_Db_class{trackerType,j} = [trueStateVTR2_All_Db_class{trackerType,j} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                            estTracks2_All_Db_class{trackerType,j}    = [estTracks2_All_Db_class{trackerType,j}    estTracks2_All{vtrDbNum,trackerType,j,class}];                

                            % reshape so that we collapse tracks across all utterances, all phonetic classes, and all formants, preserving
                            % tracker type
                            trueStateVTR2_All_Db_class_f{trackerType} = [trueStateVTR2_All_Db_class_f{trackerType} trueStateVTR2_All{vtrDbNum,trackerType,j,class}];
                            estTracks2_All_Db_class_f{trackerType}    = [estTracks2_All_Db_class_f{trackerType}    estTracks2_All{vtrDbNum,trackerType,j,class}];
                        end
                    end
                end
            end
        end

        %% Run this cell after loading VTR_error.mat, VTR_errormale.mat, VTR_errorfemale.mat or running one of previous cells 
        % Number of frames for each phonetic class
        for ii = 1:6, temp = trueStateVTR2_All_Db{1,1,ii};disp(length(temp)), end

        % Number of frames for each formant (should be the same)
        for ii = 1:3, temp = trueStateVTR2_All_Db_class{1,ii};disp(length(temp)), end

        % Number of frames for each tracker
        for ii = 1:3, temp = trueStateVTR2_All_Db_class_f{ii};disp(length(temp)), end

        % now, find error for each formant trajectory and put in table a la Deng (2006)
        % Row 1 (Class 1): F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
        % Row 2 (Class 2): F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
        % ...
        rmse_table = zeros(numClasses, numTrackers*numFormants_error);
        for class = 1:numClasses
            for trackerType = trackersToRun
                for j = 1:numFormants_error
                    differ = estTracks2_All_Db{trackerType,j,class} - trueStateVTR2_All_Db{trackerType,j,class};
                    differ = differ(~isnan(differ)); % get rid of NaNs
                    % calculate RMSE
                    rmse_table(class, j+(trackerType-1)*numTrackers) = sqrt(mean((differ).^2));

                    % calculate mean absolute error
                    %rmse_table(class, j+(trackerType-1)*numTrackers) = mean(abs(differ));
                end
            end
        end
        round(rmse_table)

        % now, find error for each formant trajectory and put in table
        % F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
        rmse_table2 = zeros(1, numTrackers*numFormants_error);
        for trackerType = trackersToRun
            for j = 1:numFormants_error
                differ = estTracks2_All_Db_class{trackerType,j} - trueStateVTR2_All_Db_class{trackerType,j};
                differ = differ(~isnan(differ)); % get rid of NaNs
                % calculate RMSE
                rmse_table2(j+(trackerType-1)*numTrackers) = sqrt(mean((differ).^2));

                % calculate mean absolute error
                %rmse_table2(j+(trackerType-1)*numTrackers) = mean(abs(differ));
            end
        end
        round(rmse_table2)

        % now, find overall error for each tracker and put in table
        % F1, F2, F3 (for all trackers)
        rmse_table3 = zeros(1, numTrackers);
        for trackerType = trackersToRun
            differ = estTracks2_All_Db_class_f{trackerType} - trueStateVTR2_All_Db_class_f{trackerType};
            differ = differ(~isnan(differ)); % get rid of NaNs
            % calculate RMSE
            rmse_table3(trackerType) = sqrt(mean((differ).^2));

            % calculate mean absolute error
            %rmse_table3(trackerType) = mean(abs(differ));
        end
        round(rmse_table3)

        % pick appropriate save function
        % save(fullfile(savedir, [prefix, '_error']))
        % save(fullfile(savedir, [prefix, '_errormale']))
        save(fullfile(savedir, [prefix, '_errorfemale']))