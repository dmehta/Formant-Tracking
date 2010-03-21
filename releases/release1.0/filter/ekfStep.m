% Do one full iteration of the Extended Kalman Filter
% Also compute likelihood of a particle
% INPUT
%
%   y - observation at iteration
%   m_up - Output of last iteration
%   P_up - Output of last iteration
%   regVar - current guess at the regularization variance
%   cepOrder - Cepstral order
%   stateVar - variance of the state
%   fs     - sampling rate
function [m_upNew, P_upNew, weight] = ekfStep(y, curInds, curBwVec, m_up, P_up, regVar, stateVar, cepOrder, fs,F_k)

% Define necessary matrices
R_k = regVar*eye(cepOrder);

m_pred = F_k*m_up;                      % Update m_{k|k-1}
P_pred = F_k*P_up*F_k' + stateVar;      % Update P_{k|k-1}

%Linearize about state
[H_k h_k] = get_linear_obserII(m_pred, curBwVec, cepOrder, fs);

S_k  = H_k*P_pred*H_k' + R_k;
y_pred = H_k*m_pred+h_k;                % Update y_{k|k-1}

% Compute particle weight
weight = my_mvnpdf(y_pred, y, S_k) + 1e-99;

% Repeat mask for observation matrix
mask = repmat(curInds,cepOrder,1);

% Update steps
gain = P_pred*H_k'*inv(S_k);              % Calculate Gain
gain = gain.*mask';                       % Mask Gain
m_upNew = m_pred + gain*(y-y_pred);       % Calculate m_{k|k}
P_upNew = P_pred - gain*H_k*P_pred;       % Calculate P_{k|k}