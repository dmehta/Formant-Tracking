% calculate exact Cramer-Rao lower bound on variance of AR coefficient estimators
% 
% see Eq. (19) [1]
% [1] B. Friedlander and B. Porat, "The exact Cramer-Rao bound for Gaussian
% autoregressive processes," IEEE Transactions on Aerospace and Electronic
% Systems, vol. 25, no. 1, pp. 3-7, 1989. 
function [CRLB, CRLB_sig, CRLB_mat, J, Rn] = AR_CRLB_exact(coefs, N, sigma2)

n = length(coefs);

% J = Jbar + (1 - n/N)*(Fisher information in asymptotic case = Jasy), Eq. (19)

% Fisher information for asymptotic case
[~, ~, ~, Jasy, Rn] = AR_CRLB2(coefs, N, sigma2);

% now all we need is Jbar

%% Although Sn = inv(Rn), solve for Sn symbolically to enable
% partial differentiation later
clear A1 A2

% define symbolic variables
for ii = 1:n
    eval(['syms a', num2str(ii), ' real'])
end

for ii = 1:n
    for jj = 1:n
        if ii == jj
            A1(ii, jj) = sym(1);
        else
            if ii > jj
                eval(['A1(ii, jj) = a', num2str(ii-jj)]);
            else
                A1(ii, jj) = sym(0);
            end
        end
    end
end

for ii = 1:n
    for jj = 1:n
        if ii >= jj
            eval(['A2(ii, jj) = a', num2str(n-ii+jj)]);
        else
            A2(ii, jj) = sym(0);
        end
    end
end

Sn = 1/sigma2 * (A1*A1' - A2*A2'); % Eq. (16)

%% solve for Jbar
%% ******************************** START HERE BY CODING DIFFERENTIATION
%% *********
Jbar = zeros(n+1,n+1);

for kk = 1:n
    for ll = 1:n
        if kk == 1 && ll == 1
            Jbar(1,1) = n/(2*sigma2^2); % Eq. (20A)
        else
            if kk == 1
                %Jbar(1,ll+1) = -1/(2*sigma2)*trace(dSn/dall*Rn); % Eq. (20B)
            else
                if ll == 1
                    %Jbar(kk+1,1) = -1/(2*sigma2)*trace(dSn/dakk*Rn); % Eq. (20B)
                else
                    %Jbar(kk+1,ll+1) = 1/2*trace(dSn/dakk*Rn*dSn/dall*Rn); % Eq. (20C)
                end
            end
        end
    end
end

%% Fisher information matrix J
J = Jbar + (1 - n/N)*Jasy;

%% CRLB of parameters {sigma2, a1, a2, ..., an}
CRLB_mat = inv(J);

CRLB_sig = CRLB_mat(1,1);
CRLB = diag(CRLB_mat(2:end, 2:end))';