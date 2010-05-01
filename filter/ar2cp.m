function cepCoeffs = ar2cp(arCoeffs, cepOrder)

p = length(arCoeffs);  % Number AR coefficients given
C = zeros(cepOrder,1); % Number Cepstral coefficients to produce
C(1) = arCoeffs(1);    % Initial case in recursion

% Case 1: index lower than max num AR coeffs
% for j = 2:p
%     a_seq = arCoeffs(j:-1:1);
%     aaInds = 1:1:j-1;
%     c_seq = C(aaInds);
%     C(j) = a_seq(1) + sum(aaInds'.*c_seq(aaInds).*a_seq(aaInds+1))/j;
% end

% Case 2: index higher than max num AR coeffs
% for j = p+1:cepOrder
%     a_seq = arCoeffs(p:-1:1)';
%     c_seq = C(j-p:1:j-1);
%     aaInds = 1:1:length(c_seq);
%     C(j) = sum(aaInds.*c_seq(aaInds)'.*a_seq(aaInds))/p;
% end

% Dan corrected 4/26/10

if cepOrder > 1
    if cepOrder < p
        jmax = cepOrder;
    else
        jmax = p;
    end
    
    for j = 2:jmax

        tmp = arCoeffs(j);
        for i = 1:j-1
            curA = arCoeffs(j-i);
            tmp = tmp + i*curA*C(i)/j;
        end
        C(j) = tmp;
    end

    if cepOrder > p
        for j = p+1:cepOrder

            tmp = 0;
            for i = j-p:j-1
                curA = arCoeffs(j-i);
                tmp = tmp + i*curA*C(i)/j;
            end
            C(j) = tmp;
        end
    end
end

cepCoeffs = C; % Set output