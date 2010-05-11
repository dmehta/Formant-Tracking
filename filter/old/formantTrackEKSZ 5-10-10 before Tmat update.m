function [m_upS, P_upS, PP_upS, lgLkl] = formantTrackEKSZ(y, F, Q, R, x0, formantInds, fs, bwStates, numF, smooth)

% This function implements a formant tracking algorithm based on an
% extended Kalman filter that tracks both poles and zeros of the vocal
% tract transfer function.
%
% Most importantly, the assumption is that the state is ordered as follows:
% x_k = [formantFrequencies; formantBandwidths; antiformantFrequencies; antiformantBandwidths]
% The bandwidths are tracked if bwStates is provided.
%
% INPUT
%   y           - Observations (complex cepstrum)
%   F           - Transition matrix
%   Q           - Process noise covariance matrix
%   R           - Observation noise covariance matrix
%   x0          - Initial State in format: [formantFrequencies;
%                   formantBandwidths; antiformantFrequencies;
%                   antiformantBandwidths]
%   formantInds - Indices for when to 'coast' formants or when not to
%   fs          - Sampling frequency, Hz
%   bwStates    - Bandwidth values. If [], then bandwidths will be tracked.
%                   Otherwise must provide bandwidth values (numFormants x
%                   numObs)
%   numF - number of formants in state, the number of antiformants and
%                   bandwidths is derived
%
% OUTPUT
%   m_upS  - Posterior means
%   P_upS  - Posterior covariances
%   PP_upS - Posterior lag-one covariances
%   lgLkl  - Log likelihood of parameters (F,Q,R,x0) given the data
%
%   Author: Daniel Rudoy, Daryush Mehta
%   Created: 02/28/2008
%   Modified: 02/10/2010, 03/21/2010
%             05/9/2010 (bandwidth tracking)

% If bandwidth tracks not provided, have to track them
trackBW = isempty(bwStates);
if trackBW
    numAntiF = size(x0,1)/2-numF;
else
    numAntiF = size(x0,1)-numF;
end

absSort = 1; %
[cepOrder N] = size(y);   % Number of cepstral coefficients and observations
lgLkl   = 0;

%Initial state and variance estimate
m_pred = F*x0;
P_pred = F*Q*F' + Q;

for k = 1:N
    
    % Linearize about state (uses Taylor expansion)
    if(trackBW)
        curFVals = m_pred([1:numF, 2*numF+1:2*numF+numAntiF],k);
        curBVals = m_pred([numF+1:2*numF, 2*numF+numAntiF+1:end],k);
        H = getH_FZBW(curFVals, curBVals, numF, cepOrder, fs);
    else
        curFVals = m_pred(:,k);
        curBVals = bwStates(:,k);
        H = getH_FZ(curFVals, curBVals, numF, cepOrder, fs);
    end
    
    curInds = formantInds(k,:); % Pull out coasting indices
    if(trackBW)
        mask = diag([curInds curInds]); % Repeat mask for observation matrix
    else
        mask = diag(curInds);
    end

    % Compute gain and knock out unobservables
    S = H*P_pred(:,:,k)*H' + R;
    K = mask*P_pred(:,:,k)*H'*inv(S); % K is the gain
     
    % Update steps
    ypred = fb2cp(curFVals(1:numF,k), curBVals(1:numF,k),cepOrder,fs)'-...
        fb2cp(curFVals(numF+1:end,k), curBVals(numF+1:end,k),cepOrder,fs)';
    
    % Calculate innovation and likelihood
    e     = y(:,k) - ypred;  % Innovation             
    lgLkl = lgLkl + gaussian_prob(e, zeros(1,length(e)), S, 1); % Likelihood
    
    % Update steps
    m_up(:,k)   = abs(m_pred(:,k) + K*e);
    P_up(:,:,k) = P_pred(:,:,k) - K*H*P_pred(:,:,k);

    % These calculations are needed for subsequent implementation of EM
    % loops
    if(k > 1)
        PP_up(:,:,k) = (eye(length(m_up(:,k)))-K*H)*F*P_up(:,:,k-1);
    else
        PP_up(:,:,k) = (eye(length(m_up(:,k)))-K*H)*F*Q;
    end
    
    % This has some regularizing effect on behaviors that are seen
    % in many synthetic examples (but not necessarily real speech examples)
    % This of course means loss of MMSE, but could potentially be replaced
    % by constraining the state with inequality constraints. Specifically
    % want to guard against going below 0Hz (or say 60 Hz) and above fs/2  
    if(absSort)
        % Sorting of the formants and zeros happens independently
        [m_up(1:numF,k) I1] = sort(abs(m_up(1:numF,k))); % formant frequencies
        if trackBW
            [m_up(2*numF+1:2*numF+numAntiF,k) I2] = sort(abs(m_up(2*numF+1:2*numF+numAntiF,k))); % anti-formant frequencies
            
            m_up(numF+1:end,k) = m_up(numF+I,k);
            
            m_up(numF+1:2*numF,k) = m_up(numF+I1,k); % formant bandwidths
            m_up(2*numF+numAntiF+1:end,k) = m_up(2*numF+numAntiF+I2,k); % anti-formant bandwidths
        else
            [m_up(numF+1:end,k) I1] = sort(abs(m_up(numF+1:end,k))); % anti-formant frequencies
        end
       
        % Build permutation matrix to adjust covariance
        Tmat = zeros(size(F));
        for i = 1:length(F)
            if(i <= numF)
                Tmat(i,I1(i)) = 1;
            else
                Tmat(i,I2(i-numF)) = 1;
            end
        end
        % Adjust covariance matrix
        P_up(:,:,k) = Tmat*P_up(:,:,k)*Tmat';
    end
  
   % Compute prediction steps for all but the last step
    if k < N
        m_pred(:,k+1)   = F*m_up(:,k);
        P_pred(:,:,k+1) = F*P_up(:,:,k)*F' + Q;
        % Below another fix, for reflection about the Nyquist frequency
        inds = m_pred(:,k+1) > fs/2;
        if(sum(inds > 0))
           m_pred(inds,k+1) = fs/2-(m_pred(inds,k+1)-fs/2);
        end
    end
    
end

% Kalman Smoothing: These are the Rauch, Tung and Striebel recursions
if(smooth)
    %These are the Rauch, Tung and Striebel recursions
    m_upS(:,N)     = m_up(:,N);
    P_upS(:,:,N)   = P_up(:,:,N);
    
    for k = (N-1):-1:1
        % Compute prediction steps for all but the last step
        sgain = P_up(:,:,k)*F'*inv(F*P_up(:,:,k)*F' + Q);
        m_upS(:,k)     = m_up(:,k)  + sgain*(m_upS(:,k+1)  - m_pred(:,k+1));
        P_upS(:,:,k)   = P_up(:,:,k)+ sgain*(P_upS(:,:,k+1) - P_pred(:,:,k+1))*sgain';
        
        % Need this calculation for subsequent EM iterations
        PP_upS(:,:,k) = PP_up(:,:,k) + (P_upS(:,:,k)-P_up(:,:,k))*inv(P_up(:,:,k))*PP_up(:,:,k);
    end
else
    % Just package the filtered outputs
    m_upS = m_up;
    P_upS = P_up;
    PP_upS = []; 
end