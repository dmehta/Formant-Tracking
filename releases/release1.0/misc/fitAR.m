% Fits AR model of order modelOrder to the signal x
% Returns mean squared error and residual
function [mse res coeffs] = fitAR(x, modelOrder, plotVar)

% Number of samples
N = length(x);

for p = modelOrder

    a = lpc(x, p);
    res = filter([0 -a(2:end)],1,x);
    for i = p+1:N
        res2(i) = 0;
        for j = 1:p
            res2(i) = res2(i) - a(j+1)*x(i-j);
        end
    end

    e = x(p+1:end) - res(p+1:end);

    %figure;
    %plot(x(p-1:end-1));
    %hold on;
    %plot(res(p:end),'g');

    % Compute MSE
    mse = mean(e.*e);
    % Compute AIC penalty
    aic = 2*p + N*log(sum(e.'*e)/N);
    % [i jbtest(e) kstest(e)]
end

if(plotVar)
    figure;
    subplot(2,1,1)
    plot(1:length(aic), aic,'r');
    xlabel('aic');
    subplot(2,1,2)
    plot(1:length(mse), mse,'b');
    xlabel('mse');
end

coeffs = -a(2:end);
