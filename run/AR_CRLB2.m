% calculate asymptotic (N -> inf) Cramer-Rao lower bound on variance of AR coefficient estimators
% 
% see Eq. (5) [1]
% [1] B. Friedlander and B. Porat, "The exact Cramer-Rao bound for Gaussian
% autoregressive processes," IEEE Transactions on Aerospace and Electronic
% Systems, vol. 25, no. 1, pp. 3-7, 1989. 
function [CRLB, CRLB_sig, CRLB_mat, J, Rn] = AR_CRLB2(coefs, N, sigma2)

n = length(coefs);

% now with covariance structure, Eq. (5) [1]

% Eq. (17):
A1 = zeros(n, n);
for ii = 1:n
    for jj = 1:n
        if ii == jj
            A1(ii, jj) = 1;
        end
        if ii > jj
            A1(ii, jj) = coefs(ii-jj);
        end
    end
end

A2 = zeros(n, n);
for ii = 1:n
    for jj = 1:n
        if ii >= jj
            A2(ii, jj) = coefs(n-ii+jj);
        end
    end
end

Rn = sigma2 .* inv(A1*A1' - A2*A2');

J = zeros(n+1,n+1);
J(1,1) = 1/(2*sigma2);
J(2:end, 2:end) = Rn;
J = N/sigma2 * J;

CRLB_mat = inv(J);

CRLB_sig = CRLB_mat(1,1);
CRLB = diag(CRLB_mat(2:end, 2:end))';