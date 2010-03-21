function [G g] = genLinearF(freq_vector, BW_vector, M, fs)

cepVec = 1:1:M;    % Index of cepstral coefficients
numFRegions = 1000; % Number of regions for formant/cepstrum linearization

regFBounds = linspace(0,2*pi, numFRegions);

% Calculations
for freqInd = 1:length(freq_vector)

    % Assign variables from state vector
    freq = freq_vector(freqInd);
    BW = BW_vector(freqInd);
    
    if(freq > fs)
        f_l = regFBounds(end-1)*fs/(2*pi);
        f_h = regFBounds(end)*fs/(2*pi);
    elseif (freq < 0)
        f_l = regFBounds(1)*fs/(2*pi);
        f_h = regFBounds(2)*fs/(2*pi);
    else
        [N f_inds] = histc(2*pi*freq/fs,regFBounds);
        f_l = regFBounds(f_inds)*fs/(2*pi);
        f_h = regFBounds(f_inds+1)*fs/(2*pi);
    end

    f_l = repmat(f_l, 1, M);
    f_h = repmat(f_h, 1, M);

        % Calculate Cepstral Coefficients
    B = 2*exp(-pi*cepVec.*BW./fs)./cepVec;
    
    y_l = B.*cos(2*pi*cepVec.*f_l/fs);
    y_h = B.*cos(2*pi*cepVec.*f_h/fs);
    
    % Slope and y-intercept
    alpha(:,freqInd) = (y_h - y_l)./(f_h - f_l);
    beta(:,freqInd) = y_l - alpha(:,freqInd)'.*f_l;  
    
end

G = alpha;
g = sum(beta,2);

