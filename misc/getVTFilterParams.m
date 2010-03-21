% Select parameters for Vocal Tract Filter
% Following Stevents p. 288
function [F, BW] = getVTFilterParams(vowel, gender)

% Female Vocal Tract Resonances
if ( gender == 'f' )
    switch vowel
        case 'i'
            F = [310 2790 3310]; % /i/ formant frequencies in Hz
            BW= [75 50 95];
        case 'e'
            F = [560 2320 2950]; % /e/ formant frequencies in Hz
            BW= [48 50 150];
        case 'ae'
            F = [860 2050 2850]; % /ae/ formant frequencies in Hz
            BW= [51 65 105];
        case 'a'
            F = [850 1220 2810]; % /a/ formant frequencies in Hz
            BW= [51 55 105];
        case 'o'
            F = [600 1200 2540]; % /o/ formant frequencies in Hz
            BW= [47 50 55];
        case 'u'
            F = [370 950 2670]; % /u/ formant frequencies in Hz
            BW= [80 50 125];
        otherwise
            error('Unknown vowel input.')
    end
    % Male vocal tract resonances, Stevens p.288
elseif ( gender == 'm' )
    switch vowel
        case 'i'
            F = [270 2290 3010]; % /i/ formant frequencies in Hz
            BW= [60 37 110];
        case 'e'
            F = [460 1890 2670]; % /e/ formant frequencies in Hz
            BW= [40 37 115];
        case 'ae'
            F = [660 1720 2410]; % /ae/ formant frequencies in Hz
            BW= [40 60  120];
        case 'a'
            F = [730 1090 2440]; % /a/ formant frequencies in Hz
            BW= [60 42 120];
        case 'o'
            F = [450 1050 2610]; % /o/ formant frequencies in Hz
            BW= [45 37 65];
        case 'u'
            F = [300 870 2240]; % /u/ formant frequencies in Hz
            BW= [60 37 50];
        otherwise
            error('Unknown vowel input.');
    end
else
    error('Gender incorrectly specified.')
end