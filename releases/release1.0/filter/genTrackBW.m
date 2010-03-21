%
% Author: Daniel Rudoy
%
% Created:  03/13/2007
% Modified: 12/13/2007

function[BW_vector] = genTrackBW(trBW_flag, BW_data)
%function[BW_vector] = genTrackBW(trBW_flag, BW_data)
% Averages, or does not average, track bandwidth information

[rows columns] = size(BW_data);

if trBW_flag == 1
    % Average Bandwidths
    BW_vector = repmat(mean(BW_data,2),1,columns);
else
    BW_vector = BW_data;
end
