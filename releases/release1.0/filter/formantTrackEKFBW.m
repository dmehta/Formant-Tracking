function [m_up, P_up] = formantTrackEKFBW(data, F, Q, R, x0, formantInds, fs)
%function [m_up, P_up] = formantTrackEKF(data, formantInds, pNoiseVar,oNoiseVar,fs,trBW_flag,BW_data,initial_state)
% Extended Kalman Filter for formant tracking (Deng et al., 2007) 
% (with indicator variables)

[cepOrder N] = size(data);
numFormants = length(x0)/2;

%Initial State and Variance Estimate
m_pred = F*x0;

P_pred = F*Q*F' + Q;

for k = 1:N    
    %Linearize about state
    [H h] = genLinearFBW(m_pred(:,k), cepOrder, fs); 

    % Pull out coasting indices
    curInds = formantInds(k,:);
    % Repeat mask for observation matrix
    mask = repmat(curInds,cepOrder,1);
    
    mask = repmat(mask,1,2);
    
    %Compute Gain 
    gain = P_pred(:,:,k)*H'*inv(H*P_pred(:,:,k)*H' + R);
    gain = gain.*mask'; % Knock out the gain as in Schmidt-Kalman filter
    
    % Update steps
    m_up(:,k) = abs(m_pred(:,k) + gain*(data(:,k)-(H*m_pred(:,k)+ h)));
    [m_up(1:numFormants,k), I] = sort(m_up(1:numFormants,k));
    P_up(:,:,k) = P_pred(:,:,k) - gain*H*P_pred(:,:,k);
    
    Tmat = zeros(size(F));
    for i = 1:numFormants
        Tmat(1,I(1)) = 1;
    end

    P_up(:,:,k) = Tmat*P_up(:,:,k)*Tmat';
    
    % Compute prediction steps for all but the last step
    if k < N
        F_k = F;
        for l = 1:length(curInds)
            if(curInds(l) == 0)
                F_k(l,:) = zeros(1,2*length(curInds));
                F_k(l+numFormants,:) = zeros(1,2*length(curInds));
                F_k(:,l) = zeros(2*length(curInds),1);
                F_k(:,l+numFormants) = zeros(2*length(curInds),1);
                F_k(l,l) = 1;
                F_k(l+numFormants,l+numFormants) = 1;
            end
        end
        m_pred(:,k+1) = F_k*m_up(:,k);
        P_pred(:,:,k+1) = F_k*P_up(:,:,k)*F_k' + Q;
    end
end
