% Input
% curState - current state vector (formants, bandwidths)
% M        - number of cepstral coefficients observed
% fs       - sampling frequency
%
%
% Created:      01/10/2008
% Last Updated: 01/18/2008
% Daniel Rudoy

function [H h] = genLinearFBW(curState, M, fs)

numFormants = length(curState)/2;
curFreq = curState(1:numFormants);
curBW   = curState(numFormants+1:end);

cepVec = 1:1:M;    % Index of cepstral coefficients

numFRegions = 1000; % Number of regions for formant/cepstrum linearization
numBRegions = 500;   % Number of regions for bandwidth/cepstrum linearization

regFBounds = linspace(0,2*pi, numFRegions);
regBBounds = linspace(0,500, numBRegions);

% Calculations
for formInd = 1:numFormants

    % Assign variables from state vector
    freq = curFreq(formInd);
    BW = curBW(formInd);
    
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
    y_l = cos(2*pi*cepVec.*f_l/fs);
    y_h = cos(2*pi*cepVec.*f_h/fs);

    % Slope and y-intercept
    alpha(:,formInd) = (y_h - y_l)./(f_h - f_l);
    beta(:,formInd) = y_l - alpha(:,formInd)'.*f_l;

    % Now linearize the bandwidth function (simultaneous over all i)
    [N b_inds] = histc(BW, regBBounds);

    if(BW < 0)
        b_l = 1;
        b_h = 2;
    elseif(BW > 500)
        b_l = regBBounds(end-1);
        b_h = regBBounds(end);
    else
        b_l = regBBounds(b_inds);
        b_h = regBBounds(b_inds+1);
    end

    z_l = exp(-pi*cepVec.*b_l/fs);
    z_h = exp(-pi*cepVec.*b_h/fs);

    f0(:,formInd) = f_l;
    b0(:,formInd) = b_l;

    gamma(:,formInd) = (z_h-z_l)./(b_h - b_l);
    delta(:,formInd) = z_l - gamma(:,formInd)'.*b_l;
end

b0m = repmat(b0, M,1);

phi = alpha.*gamma.*b0m + alpha.*delta;
psi = alpha.*gamma.*f0 + beta.*gamma;
d   = beta.*delta - alpha.*gamma.*f0.*b0m;


col  = 1./(1:1:M);
mult = repmat(col, numFormants, 1)';

% Calculate H matrix:
%H = zeros(cepOrder, 2*numFormants);
H = 2*[phi.*mult, psi.*mult];
h = 2*(col'.*sum(d,2));


