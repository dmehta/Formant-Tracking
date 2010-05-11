function [m_upS, P_upS, PP_upS, lgLkl] = formantTrackEKSZ(y, F, Q, R, x0, formantInds, fs, bwStates, numF, smooth)

% This function implements a formant tracking algorithm based on an
% extended Kalman filter that tracks both poles and zeros of the vocal
% tract transfer function
%
% Most importantly, the assumption is that the state is ordered as follows:
% x_k = [poleFrequencies zeroFrequencies]
% The bandwidths are not tracked, so bwStates provides them and that is why
% the bandwidths are not in the state vector
%
% INPUT
%   y           - Obsevations (complex cepstrum)
%   F           - Transition matrix
%   Q           - Process noise covariance matrix
%   R           - Observation noise covariance matrix
%   x0          - Initial State
%   formantInds - Mask for formant indices, to coast etc..
%   fs          - Sampling frequency
%   bwStates    - Bandwidth values
%   numF - number of formants in state, the rest are for zeros
%
% OUTPUT
%
% m_upS  - Posterior means
% P_upS  - Posterior covariances
% PP_upS - Posterior lag-one covariances
% lgLkl  - Log likelihood of parameters (F,Q,R,x0) given the data
%
%   Author: Daniel Rudoy
%   Created: 02/28/2008
%   Modified: 02/10/2010, 03/21/2010

absSort = 0;
[cepOrder N] = size(y);
lgLkl   = 0;

%Initial State and Variance Estimate
m_pred = F*x0;
P_pred = F*Q*F' + Q;

for k = 1:N
    
    % Linearize about state (uses Taylor expansion)
    H  = getH_FZ(m_pred(:,k), bwStates(:,k), numF, cepOrder, fs);
    
    curInds = formantInds(k,:); % Pull out coasting indices
    mask = diag(curInds);       % Repeat mask for observation matrix  

    % Compute gain and knock out unobservables
    S = H*P_pred(:,:,k)*H' + R;
    K = mask*P_pred(:,:,k)*H'*inv(S); % K is the gain

     
    % Update steps
    ypred = fb2cp(m_pred(1:numF,k), bwStates(1:numF,k),cepOrder,fs)'-...
        fb2cp(m_pred(numF+1:end,k), bwStates(numF+1:end,k),cepOrder,fs)';
    
    % Calculate innovation and likelihood
    e     = y(:,k) - ypred;  % Innovation             
    lgLkl = lgLkl + gaussian_prob(e, zeros(1,length(e)), S, 1); 
    
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
        [m_upF_sorted I1] = sort(abs(m_up(1:numF,k)));
        [m_upZ_sorted I2] = sort(abs(m_up(numF+1:end,k)));

        m_up(:,k) = [m_upF_sorted; m_upZ_sorted];
        % Build permutation matrix to adjust covariance
        Tmat = zeros(size(F));
        for i = 1:length(F)
            if(i <=numF)
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