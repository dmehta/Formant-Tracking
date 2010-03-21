function [G g] = get_linear_obserII(freq_vector, BW_vector, cep_order, fs)
% Provide Cepstral estimate
% Given BANDWIDTH and prediction of frequency, computes the LPC cepstra up
% to order M, this is the function of the script linear_observ.mat

% Created : 3/1/2007
% Last Updated : 01/10/2007
% Dan Spendley, Daniel Rudoy

M = cep_order;     % Number of cepstral coefficients observed
cepVec = 1:1:M;    % Index of cepstral coefficients
numRegions = 1000; % Number of linearizing terms 

new_reg_bound = (1:1:numRegions)*(2*pi/numRegions) - (2*pi/numRegions);
reg_bound_rads = zeros(2,numRegions);
reg_bound_rads(1,:) = new_reg_bound;
reg_bound_rads(2,1:end-1) = new_reg_bound(2:end); 
reg_bound_rads(2,end) = 2*pi;

% Calculations
for formInd = 1:length(freq_vector)
    
    % Assign variables from state vector
    freq = freq_vector(formInd);
    BW = BW_vector(formInd);
        
    phase_estimate = mod((2*cepVec*pi*freq/fs),2*pi);
    unwrap_est = floor((2*cepVec*pi*freq/fs)/(2*pi));
    inds = (phase_estimate == 0);
    est_region(inds) = 1;       
    [N bins] = histc(phase_estimate(~inds), reg_bound_rads(1,:));
    est_region(~inds) = bins;
     
     % Need to be careful here!
    boundBinInd = bins == 0;
    est_region(boundBinInd) = 1000;
     
    % Cacluate Frequency Boundary Points
    f_l = (reg_bound_rads(1,est_region)+2*pi*unwrap_est)*fs./(2*pi*cepVec);
    f_h = (reg_bound_rads(2,est_region)+2*pi*unwrap_est)*fs./(2*pi*cepVec);

    % Calculate Cepstral Coefficients
    B = 2*exp(-pi*cepVec.*BW./fs)./cepVec;
    c_l = B.*cos(2*pi*cepVec.*f_l/fs);
    c_h = B.*cos(2*pi*cepVec.*f_h/fs);
    
    % Slope and y-intercept
    alpha(:,formInd) = (c_h - c_l)./(f_h - f_l);
    beta(:,formInd) = c_l - alpha(:,formInd)'.*f_l;  
end

G = alpha;
g = sum(beta,2);

