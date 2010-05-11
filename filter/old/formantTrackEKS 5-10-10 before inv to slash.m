function [m_upS, P_upS, PP_upS,lgLkl] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth)

% INPUT
%
% y  - Matrix of observed data (R^{num cepstral coeffs x num observations }
% F  - Transition matrix
% Q  - Process noise covariance
% R  - Observation noise covariance
% x0 - Initial state
% formantInds - Indices for when to 'coast' formants or when not to
% fs - sampling frequency (e.g., 16000)
% bwStates - If [], then bandwidths will be tracked. Otherwise must provide 
  %Bandwidth values (numFormants x num observations)
% smooth   - if 1, then run Extended Kalman Smoother, if 0 only EK Filter.

% OUTPUT
% m_upS - numF x numObservations MMSE estimate
% P_upS - numF x numF x numObservations covariances

% Author: Daniel Rudoy
% Modified: 02/14/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If bandwidth tracks not provided, have to track them
trackBW = isempty(bwStates);

absSort = 1;
[cepOrder N] = size(y);    % Number of cepstral coefficients and data points
lgLkl   = 0;

if(trackBW)
    numF  = length(x0)/2; % Number of formants (half the params)
else
    numF  = length(x0);   % Number of formants (all the params)
end

%Initial state and variance estimate
m_pred = F*x0;        % Compute m_{1,0}
P_pred = F*Q*F' + Q;  % Compute P_{1,0}

for k = 1:N

    % Linearize about m_{k,k-1} using Taylor expansion (not piecewise lin.)
    % Deng's linearization also implemented (see getLinearH)
    if(trackBW)
        curFVals = m_pred(1:numF,k);
        curBVals = m_pred(numF+1:end,k);
        H = getH_FBW(curFVals, curBVals, numF, cepOrder, fs);
    else
        curFVals = m_pred(:,k);
        curBVals = bwStates(:,k);
        H = getH_F(curFVals, curBVals, numF, cepOrder, fs);
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
    ypred = fb2cp(curFVals, curBVals, cepOrder, fs)';
    
    % Calculate innovation and likelihood
    e     = y(:,k) - ypred;               
    lgLkl = lgLkl + gaussian_prob(e, zeros(1,length(e)), S, 1); 
   
    % Update steps
    m_up(:,k)   = m_pred(:,k) + K*e;
    P_up(:,:,k) = (eye(length(m_up(:,k)))-K*H)*P_pred(:,:,k);

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
        % Lower bound and label switching fixes
        [m_up(1:numF,k), I] = sort(abs(m_up(1:numF,k)));
        Tmat = zeros(size(F));
        for i = 1:numF
            Tmat(i,I(i)) = 1;
            if(trackBW)
                Tmat(i+numF, I(i)+numF) = 1;
            end
        end
        P_up(:,:,k) = Tmat*P_up(:,:,k)*Tmat';
    end

    % Compute prediction steps for all but the last step
    if k < N
        m_pred(:,k+1)   = F*m_up(:,k);
        P_pred(:,:,k+1) = F*P_up(:,:,k)*F' + Q;
        
        inds = m_pred(:,k+1) > fs/2;
        if(sum(inds > 0))
           m_pred(inds,k+1) = fs/2-(m_pred(inds,k+1)-fs/2);
        end
        
    end
end

% Run Kalman Smoother
if(smooth)
    %These are the Rauch, Tung and Striebel recursions
    m_upS(:,N)     = m_up(:,N);
    P_upS(:,:,N)   = P_up(:,:,N);
    
    for k = (N-1):-1:1
        % Compute prediction steps for all but the last step
        sgain = P_up(:,:,k)*F'*inv(F*P_up(:,:,k)*F' + Q);
        m_upS(:,k)     = m_up(:,k)  + sgain*(m_upS(:,k+1)  - m_pred(:,k+1));
        P_upS(:,:,k)   = P_up(:,:,k)+ sgain*(P_upS(:,:,k+1) - P_pred(:,:,k+1))*sgain';
        
        % This I think is what the Kevin Murphy version looks like sane as Ostendorf paper
        PP_upS(:,:,k) = PP_up(:,:,k) + (P_upS(:,:,k)-P_up(:,:,k))*inv(P_up(:,:,k))*PP_up(:,:,k);
        
    end
else
    % Just package the filtered outputs
    m_upS = m_up;
    P_upS = P_up;
    PP_upS = []; 
end