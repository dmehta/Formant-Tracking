function [m_upS, P_upS] = formantTrackEKSZ(y, F, Q, R, x0, formantInds, fs, bwStates, numF, smooth)

% This function implements a formant tracking algorithm based on an
% extended Kalman filter that tracks both poles and zeros of the vocal
% tract transfer function
%
% Input
%   y        - Obsevations (complex cepstrum)
%   F           - Transition matrix
%   Q           - Process noise covariance matrix
%   R           - Observation noise covariance matrix
%   x0          - Initial State
%   formantInds - Mask for formant indices, to coast etc..
%   fs          - Sampling frequency
%   bwStates    - Bandwidth values
%   numF - number of formants in state, the rest are for zeros
%
%   Author: Daniel Rudoy
%   Created: 02/28/2008
%   Modified: 02/10/2010

absSort = 0;
[cepOrder N] = size(y);

%Initial State and Variance Estimate
m_pred = F*x0;
P_pred = F*Q*F' + Q;

for k = 1:N
    %Linearize about state
    H  = getH_FZ(m_pred(:,k), bwStates(:,k), numF, cepOrder, fs);
    curInds = formantInds(k,:); % Pull out coasting indices
    mask = diag(curInds);       % Repeat mask for observation matrix  
    gain = mask*P_pred(:,:,k)*H'*inv(H*P_pred(:,:,k)*H' + R); % Compute Gain

    % Update steps
    ypred = fb2cp(m_pred(1:numF,k), bwStates(1:numF,k),cepOrder,fs)'-...
        fb2cp(m_pred(numF+1:end,k), bwStates(numF+1:end,k),cepOrder,fs)';
    
    m_up(:,k)   = m_pred(:,k) + gain*(y(:,k)-ypred);
    P_up(:,:,k) = P_pred(:,:,k) - gain*H*P_pred(:,:,k);

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
    end
end

% Kalman Smoothing: These are the Rauch, Tung and Striebel recursions
if(smooth)
    m_upS(:,N)   = m_up(:,N);
    P_upS(:,:,N) = P_up(:,:,N);

    for k = (N-1):-1:1
        sgain = P_up(:,:,k)*F'*inv(F*P_up(:,:,k)*F' + Q);
        m_upS(:,k)   = m_up(:,k)   + sgain*(m_upS(:,k+1) - m_pred(:,k+1));
        P_upS(:,:,k) = P_up(:,:,k) + sgain*(P_upS(:,:,k+1)-P_pred(:,:,k+1))*sgain';
    end
else
    m_upS = m_up;
    P_upS = P_up;
end