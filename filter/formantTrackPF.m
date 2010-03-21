function [pmean, phist, Neff] = formantTrackPF(y, F, Q, R, x0, formantInds, fs, BW_data, N)

% Particle filtering approach to tracking through nonlinearity. This is
% not used anymore, because initial experiments indicated that the EKF
% does a fine job on its own tracking through the nonlinearity.

% Author: Daniel Rudoy
% Created: 02/28/2007
% Modified: 02/23/2010

% Get implicit dimensions from input arguments
[cepOrder numIter] = size(y); numF = length(x0);

proposalEKF = 0; % This is not yet complete

% Pull out EKF matrices
if(proposalEKF)
    [m_upS, P_upS, PP_upS, lgLkl] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, BW_data, 0);
end

% Set flag for reflection of particles about 0 and Nyquist
reflect = 1;

% Initialize particles
samples = repmat(x0, 1, N);
if reflect
    % Reflect about origin and Nyquist freq
    samples = abs(samples);
    inds = find(samples > fs/2);
    samples(inds) = fs/2-mod(samples(inds),fs/2);
end

% Create/initialize needed loop variables
allSamples = zeros(numIter, numF, N);
pmean = zeros(numF,numIter);

lweights = zeros(1,N);
for k = 1:numIter

    allSamples(k,:,:) = samples;     % Store samples


    % Compute likelihood of every particle
    curObs = y(:,k);
    for p = 1:N

        % Map each particle into cepstral domain
        ypred = fb2cp(samples(:,p), BW_data(:,k), cepOrder, fs)';
        
        % Compute log-likelihood of particle
        diffSq = (ypred-curObs)'*inv(R)*(ypred-curObs);
        curlklWt = -.5*sum(diffSq) -.5*log(det(R))+ 1e-99;
        
        if(proposalEKF)
            priorWt = 0; % not yet coded
            propWt  = 0; % not yet coded
            lweights(1,p) = lweights(1,p) + curlklWt + priorWt - propWt;
        else
            lweights(1,p) = lweights(1,p) + curlklWt;
        end
    end

    maxLWeight = max(lweights);           % Find max weights
    lweights = lweights + (1-maxLWeight); % Shift weights up
    weights = exp(lweights);              % Exponentiate
    weights = weights/sum(weights);       % Normalize weights

    % Compute posterior mean
    pmean(:,k) = samples*weights';

    % Check number of effective particles and resample if it is < N/2
    Neff(k) = 1/sum(weights*weights');
    if(Neff(k) < floor(N/2))
        % Resample particles need to replace with stratified sampling later
        %newInds = randsample(N,N,true,weights); % Multinomial Resampling
        newInds = systematicR(1:1:N, weights');  % Systematic Resampling
        samples(:,1:N) = samples(:,newInds);
        lweights = zeros(1,N);   % Reset log weights to 0
    else
        % Do nothing
    end

    % Evolve samples according to model
    sQ = sqrtm(Q);
    for p = 1:N
        curFVals = samples(:,p);
        curBVals = BW_data(:,k);
        
        % Evolution to according to prior
        m_pred = F*curFVals;
        if(proposalEKF & k > 2)
            %% This part is coded, but the weight update is not, because
            %% we need to update all the covariances etc...
            H = getH_F(curFVals, curBVals, numF, cepOrder, fs);
            P_pred = F*P_upS(:,:,k-1)*F' + Q;
            S = H*P_pred*H' + R;
            K = P_pred*H'*inv(S);
            
            ypred = fb2cp(curFVals, curBVals, cepOrder, fs)';
            samples(:,p) = m_pred + K*(y(:,k)-ypred);
        else
            samples(:,p) = (m_pred + sQ*randn(numF,1));
        end
    end
    if reflect
        % Reflect about origin and Nyquist freq
        samples = sort(abs(samples));
        inds = find(samples>fs/2);
        samples(inds) = fs/2-mod(samples(inds),fs/2);
    end
end
phist = []; %TODO: histogram posterior
