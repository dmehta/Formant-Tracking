function [g_matrix g_vector] = get_linear_obserII(freq_vector, BW_vector, cep_order, f_samp)
%function [g_matrix g_vector] = get_linear_obserII(freq_vector, BW_vector, cep_order, f_samp)
% Provide Cepstral estimate
% Given BANDWIDTH and prediction of frequency, computes the LPC cepstra up
% to order M, this is the function of the script linear_observ.mat

% Created : 3/1/2007
% Last Updated : 3/16/2007
% Dan Spendley, Daniel Rudoy

% Variable Declarations
M = cep_order;
num_regions = 1000;

new_reg_bound = (1:1:num_regions)*(2*pi/num_regions) - (2*pi/num_regions);
reg_bound_rads = zeros(2,num_regions);
reg_bound_rads(1,:) = new_reg_bound;
reg_bound_rads(2,1:end-1) = new_reg_bound(2:end); 
reg_bound_rads(2,end) = 2*pi;

% Calculations
for pole_index = 1:length(freq_vector)
    
    % Assign variables from state vector
    freq = freq_vector(pole_index);
    BW = BW_vector(pole_index);
    
    iiVect = 1:1:M;
    B = 2*exp(-pi*iiVect.*BW./f_samp)./iiVect;
    
    phase_estimate = mod((2*iiVect*pi*freq/f_samp),2*pi);
    unwrap_est = floor((2*iiVect*pi*freq/f_samp)/(2*pi));
    est_region2 = zeros(1,M); 
    inds = (phase_estimate == 0);
    est_region(inds) = 1;       
    [N bins] = histc(phase_estimate(~inds), reg_bound_rads(1,:));
    est_region(~inds) = bins;
     
     % Need to be careful here!
    boundBinInd = bins == 0;
    est_region(boundBinInd) = 1000;
     
    % Cacluate Frequency Boundary Points
    f_l = (reg_bound_rads(1,est_region)+2*pi*unwrap_est)*f_samp./(2*pi*iiVect);
    f_h = (reg_bound_rads(2,est_region)+2*pi*unwrap_est)*f_samp./(2*pi*iiVect);

    % Calculate Cepstral Coefficients
    c_l = B.*cos(2*pi*iiVect.*f_l/f_samp);
    c_h = B.*cos(2*pi*iiVect.*f_h/f_samp);
    
    % Slope and y-intercept
    alpha(:,pole_index) = (c_h - c_l)./(f_h - f_l);
    beta(:,pole_index) = c_l - alpha(:,pole_index)'.*f_l;  
end

% Calculate intercept terms
for jj = 1:M
    g2(jj) = sum(beta(jj,:));
end

% Output correct variables
g_matrix = alpha;
g_vector = g2';
