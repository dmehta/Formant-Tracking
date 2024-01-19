        %% For JASA re-submission: Calculate error with no class division
        %% (e.g., for VTRSynth)
        trueStateVTR2_All_Db = cell(numTrackers, numFormants_error);
        estTracks2_All_Db = cell(numTrackers, numFormants_error);
        trueStateVTR2_All_Db_f = cell(1, numTrackers);
        estTracks2_All_Db_f = cell(1, numTrackers);
        for trackerType = trackersToRun
            for j = 1:numFormants_error
                for vtrDbNum = vtrDbNums % *** set to vtrDbNum if only interested in current utterance RMSE (otherwise, vtrDbNums)
                %for vtrDbNum = [1:329 331:516] % VTR F3 blatantly incorrect so skews results    
                    % reshape so that we collapse tracks across all utterances, preserving
                    % categories of formant number and tracker type
                    trueStateVTR2_All_Db{trackerType,j} = [trueStateVTR2_All_Db{trackerType,j} trueStateVTR2_All{vtrDbNum,trackerType,j}];
                    estTracks2_All_Db{trackerType,j}    = [estTracks2_All_Db{trackerType,j}    estTracks2_All{vtrDbNum,trackerType,j}];

                    % reshape so that we collapse tracks across all utterances and all formants, preserving
                    % tracker type
                    trueStateVTR2_All_Db_f{trackerType} = [trueStateVTR2_All_Db_f{trackerType} trueStateVTR2_All{vtrDbNum,trackerType,j}];
                    estTracks2_All_Db_f{trackerType}    = [estTracks2_All_Db_f{trackerType}    estTracks2_All{vtrDbNum,trackerType,j}];                
                end
            end
        end

        % Run for no phonetic categories, after loading .mat file (e.g., VTR516.mat or VTRsynthf0516.mat) and running previous cell
        % Number of frames for each formant (should be the same)
        for ii = 1:3, temp = trueStateVTR2_All_Db{1,ii};disp(length(temp)), end

        % Number of frames for each tracker
        for ii = 1:3, temp = trueStateVTR2_All_Db_f{ii};disp(length(temp)), end

        % now, find error for each formant trajectory and put in table
        % F1, F2, F3  F1, F2, F3  F1, F2, F3 (for all trackers)
        rmse_table2 = zeros(1, numTrackers*numFormants_error);
        for trackerType = trackersToRun
            for j = 1:numFormants_error
                differ = estTracks2_All_Db{trackerType,j} - trueStateVTR2_All_Db{trackerType,j};
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
            differ = estTracks2_All_Db_f{trackerType} - trueStateVTR2_All_Db_f{trackerType};
            differ = differ(~isnan(differ)); % get rid of NaNs
            % calculate RMSE
            rmse_table3(trackerType) = sqrt(mean((differ).^2));

            % calculate mean absolute error
            %rmse_table3(trackerType) = mean(abs(differ));
        end
        round(rmse_table3)

        save(fullfile(savedir, [prefix, '_error']))
