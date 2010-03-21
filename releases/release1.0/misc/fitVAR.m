% *******************************************************
% VECTOR AR
% *******************************************************
%function [mse res PI] = fitVAR(x, p, plotVar,avg)
function [PI] = fitVAR(x, p)

% Number of points per time-series, number of time-series
[N dim] = size(x);

% Construct VAR(p) estimate following Hamilton 11.1
m1 = zeros(dim,dim*p);
m2 = zeros(dim*p);
for i = p+1:N
    yt = x(i,:).';
    xt = [];
    for j = 1:p
        xt = [xt; x(i-j,:).'];
    end
    m1 = m1 + yt * xt.';
    m2 = m2 + xt * xt.';
end

PI = m1*inv(m2);


% c = [0 0]';
% for i = 1:p
%     tmp{i} = PI(:,2*i-1:2*i);
% end

% for i = 1:p
%     tmp{i} = PI(:,3*i-1:3*i+1);
% end

% Compute Residual Errors

% for i = p+1:N
%     est(:,i) = zeros(2,1);
%     for j = 1:p
%         est(:,i) = est(:,i) + tmp{j}*x(i-j,:).';
%     end
%     err(:,i) =  x(i,:).' - est(:,i);
% end
% 
% if plotVar
%     figure;
%     subplot(2,1,1);
%     hist(err(1,:), 30);
%     axis([-1 1 0 200]);
%     subplot(2,1,2);
%     hist(err(2,:), 30);
%     axis([-1 1 0 200]);
% %     subplot(3,1,3);
% %     hist(err(3,:), 30);
% %     axis([-1 1 0 200]);
% end

%res = err'; % Store residual

% Compute MSE
% mse = [];
% for i = 1:dim
%     mse = [mse var(err(i,:))];
% end
