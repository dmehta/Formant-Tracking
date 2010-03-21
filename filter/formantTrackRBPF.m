function [xhat, errVar,pEst] = formantTrackRBPF(N, data, formantInds, pNoiseVar,oNoiseVar,fs,trBW_flag,BW_data,initial_state,F)

% Rao-Blackwell particle filter for inferring the regularization parameter
% in the Deng Model. The basic implementation will 'change' the model in
% that it will allow the regularization parameter to vary over time. Hence
% this model is slightly more powerful, and can capture stronger residual
% information. 
%
% We have not followed through on careful study of the performane in this
% model.

% Author: Daniel Rudoy

% Created: 03.25.2007

[cepOrder numIter] = size(data);
stateSize = length(initial_state);

%%%%%%%%%%%%%%% Initialization for Particle Filter Part %%%%%%%%%%%%%%%
lweights = zeros(1,N);

% Samples hold the one dimensional regularization parameter
samples = abs(.005*randn(1,N) + oNoiseVar*ones(1,N));
allSamples = zeros(numIter, N);
allWeights = zeros(numIter, N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% Extended Kalman Filter Part %%%%%%%%%%%%%%%%%%%%%%%%%%%

BW_vector = genTrackBW(trBW_flag,BW_data); % Modify BW for linearization

% Initial covariance estimates
ind_mat = eye(length(pNoiseVar));
if length(pNoiseVar) == 1
    stateVar = pNoiseVar*eye(stateSize);
else
    for ii = 1:length(pNoiseVar)
        stateVar(ii,:) = pNoiseVar(ii).*ind_mat(ii,:);
    end
end

% Initialize Kalman Matrices
m_init = initial_state;            % Initialize m_{0|0}(theta_1) = m0(theta_0)
P_init = F*stateVar*F' + stateVar; % Initialize P_{0|0}(theta_0) = P_0(theta_0)

% Initialize variables for storing means and variances
m = zeros(N, stateSize, numIter);
P = zeros(N, stateSize, stateSize*numIter);

for k = 1:numIter

    allSamples(k,:) = samples;
    curInds = formantInds(k,:);
    
    % For each sample compute the weight by doing one iteration of the EKF
    % we will use the EKF of Deng with appropriate linearization etc.. in
    % each of these steps. 
    for p = 1:N       
        regVar = samples(p);
        if(k ~= 1)
            m_up = (m(p,:,k-1)');
            tmp = P(p,:,(k-2)*stateSize+1:(k-1)*stateSize);
            P_up = reshape(tmp,stateSize,stateSize);
        else
            m_up = m_init;
            P_up = P_init;
        end
       
        [m_up P_up curWeight] = ekfStep(data(:,k), curInds, BW_vector(:,k), m_up, P_up, regVar, stateVar, cepOrder, fs,F);
        
        % Update the weights (adding because in log domain)
        if(sum(curInds) == length(curInds))
            % Only update if nothing happened with weights
            lweights(1,p) = lweights(1,p) + log(curWeight);
        else
            % Do nothing
        end
        % Store only update matrices for this particle
        m(p,:,k) = m_up';
        P(p,:,(k-1)*stateSize+1:k*stateSize) =  P_up;
    end

    %Find max weight
    maxLWeight = max(lweights);
    % Shift weights up
    lweights = lweights + (1-maxLWeight);
    % Exponentiate
    weights = exp(lweights);
    % Normalize weights
    weights = weights/sum(weights);
    allWeights(k,:) = weights;
    
    estParam(k) = samples*weights';
    
    % Check number of effective particles and resample if it is < N/2
    Neff(k) = 1/sum(weights*weights');
    if(Neff(k) < N/2)
         % Resample particles need to replace with stratified sampling later
         %newInds = randsample(N,N,true,weights); % Multinomial Resampling
         newInds = systematicR(1:1:N, weights');  % Systematic Resampling
         samples(:,1:N) = samples(:,newInds);
         weights = 1/N*ones(1,N); % Reset weights to 1/N
         lweights = zeros(1,N);   % Reset log weights to 0
     else
         % Do nothing
    end
 
     % Propagate the particles model in the log domain
     if(sum(curInds) == length(curInds))
        %lsamples = log(samples);
        %samples = exp(lsamples + sqrt(.0005)*randn(1,N));
        samples = abs(samples + sqrt(.005)*randn(1,N)); 
     end
end

% Package output
pmean = zeros(numIter, stateSize);
pmax = zeros(numIter,stateSize);
for k = 1:numIter
    for p = 1:N
        pmean(k,:) = pmean(k,:) + allWeights(k, p)*squeeze(m(p,:,k));
    end
    [y ind] = max(allWeights(k,:));
    pmax(k,:) = squeeze(m(ind,:,k));
end

xhat = pmean';     % Posterior Mean
%xhat = pmax';     % MAP
errVar = P;        % Covariances
pEst = mean(estParam);
