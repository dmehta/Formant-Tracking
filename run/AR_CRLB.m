% calculate asymptotic (N -> inf) Cramer-Rao lower bound on variance of AR coefficient estimators
% 
% see Eq. (34) [1]
% [1] B. Friedlander, "On the computation of the Cramer-Rao bound for ARMA
% parameter estimation," IEEE Transactions on Acoustics, Speech, and Signal
% Processing, vol. 32, no. 4, pp. 721-727, 1984.

function CRLB = AR_CRLB(coefs, N)

n = length(coefs);

CRLB = zeros(1, n);

if n == 1
    CRLB(1) = 1/N*(1 - coefs.^2);
    return;
end



for coef_index = 1:n
    sum1 = 1;
    sum2 = 0;
    
    if coef_index == 1
        sum1 = 1;
    else
        for jj = 1:coef_index-1
            sum1 = sum1 + coefs(jj)^2;
        end
    end

    for jj = n:-1:n-coef_index+1
        sum2 = sum2 + coefs(jj)^2;
    end
    
    CRLB(coef_index) = 1/N*(sum1 - sum2);
end