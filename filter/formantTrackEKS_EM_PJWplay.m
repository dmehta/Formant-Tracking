function [m_upS, P_upS, theta, lgLkl] = formantTrackEKS_EM(y, thetaInit, formantInds, fs, bwStates, maxNumIter)

% Extended Kalman Filter for tracking vocal tract resonances and bandwidths
% together with EM steps in order to estimate the system parameters
% theta = (F, Q, R, x0, Sigma0)

% INPUT

% OUTPUT

% Author: Daniel Rudoy
% Created:  02/16/2010
% Modified: 02/16/2010

F =  thetaInit.F;
Q =  thetaInit.Q;
R =  thetaInit.R;
x0 = thetaInit.x0;

i = 1;
while (i <= maxNumIter)

    % E-step: run Kalman smoother recursions
    smooth = 1; % Always run Kalman smoother in this setting
    [m_upS P_upS PP_upS lgLkl(i)] = formantTrackEKS(y, F, Q, R, x0, formantInds, fs, bwStates, smooth);

    % M-step: re-estimate signal parameters
    [theta] = EKS_EM(y, m_upS, P_upS, PP_upS, fs, bwStates, formantInds, thetaInit);
    % Not much of a difference, everywhere ~10Hz RMSE average
    %[theta] = EKS_EM(y, m_upS, P_upS, PP_upS, fs, bwStates, []); % Don't coast re-estimates

    F  = theta.F;
    Q  = theta.Q;
    R  = theta.R;
    %x0 = theta.x0;

    i = i + 1;
end

% figure;
% plot(lgLkl);
% title('Log Likelihood');


function [theta] = EKS_EM(y, m_upS, P_upS, PP_upS, fs, bwStates, formantInds, thetaInit)

trackBw = isempty(bwStates); % Check if bandwidths are being tracked or not
if(trackBw)
    numFormants = length(m_upS(:,1))/2;
else
    numFormants = length(m_upS(:,1));
end


[Nc, N] = size(y);  % Store number of cepstral coeffs/observations

A = zeros(size(P_upS(:,:,1))); % Initialize A
B = zeros(size(P_upS(:,:,1))); % Initialize B
C = zeros(size(P_upS(:,:,1))); % Initialize C

D = zeros(Nc,Nc);
E = zeros(Nc,Nc);
F = zeros(Nc,Nc);

% If formant indices are passed in, then do not re-estimate parameters over 
% the noise indices
useInds = isempty(formantInds); 

for t = 1:N         

    if(useInds & formantInds(t,1)==0)
        continue;
    end
    
    %% Compute all terms necessary to update F, Q
    if( t > 1 )
        A = A + P_upS(:,:,t-1) + m_upS(:,t-1)*m_upS(:,t-1)';
        B = B + PP_upS(:,:,t-1) + m_upS(:,t)*m_upS(:,t-1)'; 
        C = C + P_upS(:,:,t) + m_upS(:,t)*m_upS(:,t)';
    end   
    
    %Compute all terms necessary to update R: relinearize about the smoothed mean
    if(trackBw)
        curFVals = m_upS(1:numFormants,t); curBVals = m_upS(numFormants+1:end,t);
        H_t = getH_FBW(curFVals, curBVals, numFormants, Nc, fs);
        h_t = fb2cp(curFVals, curBVals, Nc, fs)' - H_t*m_upS(:,t);
    else
        H_t = getH_F(m_upS(:,t), bwStates(:,t), numFormants, Nc, fs);
        h_t = fb2cp(m_upS(:,t), bwStates(:,t), Nc, fs)' - H_t*m_upS(:,t);
    end
    
    yCur = y(:,t) - h_t; % Correct for the linear term

    D = D + (H_t*(P_upS(:,:,t) + m_upS(:,t)*m_upS(:,t)')*H_t');
    E = E + H_t*m_upS(:,t)*yCur';
    F = F + yCur*yCur';
end

NS = sum(formantInds(:,1)); % If nothing omitted due to silences, this 
                            % term is equal to N
A = A/(NS-1); B = B/(NS-1); C = C/(NS-1);

% Now compute different estimators of R
% Note to user: full rank R requires lots of data to get robust estimator
% A diagonal covariance seems reasonable, the most `robust' is the multiple
% of identity model
R = (F - 2*E + D)/N;               % Here R is unconstrained
Rdiag = diag(diag(R));             % Here R is diagonal                 
sigmaR = trace(R/Nc);        % Compute variance for multiple of identity model
RdiagConst = sigmaR*eye(Nc) ; % Here R is identity times a constant
RdiagDecrease = trace(R*inv(diag(1./(1:Nc)))/Nc)
RdiagDecrease = RdiagDecrease*diag(1./(1:Nc));

theta.F = B*inv(A);             % Update estimate of F
%theta.F = trace(theta.F/size(thetaInit.F,1))*eye(size(thetaInit.F));
theta.Q = (C - B*inv(A)*B');     % Update estimate of Q
theta.R = RdiagDecrease;           % Update estimate of R 
theta.x0 = m_upS(:,1);          % Update estimate of x0