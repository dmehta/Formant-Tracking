% addpath(genpath('../')); % Paths

f = [500 1000 1400];
bw = [80 120 160];
Z= []; Zbw = [];
fs = 8000;
numTrials = 100;

[N, D] = fb2tf(f, bw, Z, Zbw, fs);
 
figure, freqz(sum(D), D, 1000, fs) % set DC to 0 dB
ylim([-40 40])

% CRLBs of AR coefficients
coefs = D(2:end);
p = length(f)*2;
N = 1000;
sigma2 = 1;

CRLB = AR_CRLB(coefs, N)
CRLB2 = AR_CRLB2(coefs, N, sigma2)
CRLB_exact = AR_CRLB_exact(coefs, N, sigma2)

%%
fprime = zeros(length(f), numTrials);
bwprime = zeros(length(bw), numTrials);
aprime = zeros(2*length(f), numTrials);

randn('state',sum(100*clock)); % Seeds

for ii = 1:numTrials
    s = filter(sum(D), D, randn(1, 1000));

    a = arcov(s, length(f)*2);

    [z, p, k] = tf2zp(1, a);
    kk = find(imag(p) == 0);
    p(kk) = [];
    p = p(1:2:end);

    fp = angle(p)/(2*pi)*fs;
    bwp = -log(abs(p))*fs/pi;
    
    [fp, jj] = sort(fp);
    bwp = bwp(jj);
    
    fprime(:, ii) = fp;
    bwprime(:, ii) = bwp;
    aprime(:, ii) = a(2:end);
end

%%
figure, hold on
for ii = 1:length(f)
    disp(['Formant #', num2str(ii)])
    
    plot(fprime(ii,:), 'g.')
    plot(bwprime(ii,:),'kx')

    plot(1:numTrials, f(ii)*ones(1, numTrials), 'w')
    plot(1:numTrials, bw(ii)*ones(1, numTrials), 'w')

    disp(['Freq Bias: ', num2str(mean(fprime(ii,:))-f(ii))])
    disp(['BW Bias: ', num2str(mean(bwprime(ii,:))-bw(ii))])
    
    disp(['Freq SD: ', num2str(std(fprime(ii,:)))])
    disp(['BW SD: ', num2str(std(bwprime(ii,:)))])

    disp(['Freq Var: ', num2str(var(fprime(ii,:)))])
    disp(['BW Var: ', num2str(var(bwprime(ii,:)))])
    
    disp(['Freq RMSE: ', num2str(sqrt(mean((fprime(ii,:)-f(ii)).^2)))])
    disp(['BW RMSE: ', num2str(sqrt(mean((bwprime(ii,:)-bw(ii)).^2)))])
    
    disp(' ')
end

% ARerr = mean(aprime(ii*2-1:ii*2,:), 2)-D(ii*2:ii*2+1)';
% disp(['AR Coef Bias: ', num2str(ARerr')])
% disp(['AR Coef SD: ', num2str(std(aprime(ii*2-1:ii*2,:), 0, 2)')])
% disp(['AR Coef Var: ', num2str(var(aprime(ii*2-1:ii*2,:), 0, 2)')])
% 
% err = aprime(ii*2-1:ii*2,:)-repmat(D(ii*2:ii*2+1)', 1, numTrials);
% disp(['AR Coef RMSE: ', num2str(sqrt(mean(err.^2, 2)'))])    

ARerr = mean(aprime, 2)-D(2:end)';
disp(['AR Coef Bias: ', num2str(ARerr')])
disp(['AR Coef SD: ', num2str(std(aprime, 0, 2)')])
disp(['AR Coef Var: ', num2str(var(aprime, 0, 2)')])

err = aprime-repmat(D(2:end)', 1, numTrials);
disp(['AR Coef RMSE: ', num2str(sqrt(mean(err.^2, 2)'))])

format_plot
ylim([0 max(f)+500])

%% given CRLB as variance, propagate through equations for single resonator
clear

f = [500];
bw = [80];
Z= []; Zbw = [];
fs = 8000;
numTrials = 100;

[N, D] = fb2tf(f, bw, Z, Zbw, fs);
 
figure, freqz(sum(D), D, 1000, fs) % set DC to 0 dB
ylim([-40 40])

% CRLBs of AR coefficients
coefs = D(2:end);
p = length(f)*2;
N = 1000;

CRLB = AR_CRLB(p, N, coefs);

% AR to center frequency/bandwidth, formant equation
% start here by generating multiple sequences
% randn('state',sum(100*clock)); % Seeds
randn('state',1); % Seeds

coef1_pdf = coefs(1) + sqrt(CRLB(1))*randn(1, 1000);
coef2_pdf = coefs(2) + sqrt(CRLB(2))*randn(1, 1000);

fprime = fs/(2*pi).*acos( (-coef1_pdf) ./ (2*sqrt(abs(coef2_pdf))) );
bwprime = -fs/pi.*log( sqrt(abs(coef2_pdf)) );

%% plot histograms of formant frequencies and bandwidths

for ii = 1:length(f)
    figure

    subplot(211)
    edges = f(ii)-100:2:f(ii)+100;
    NN = histc(fprime(ii,:), edges);
    bar(edges, NN,'histc')
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','none')
%     set(gca, 'YLim', [0 1200])
%     set(gca, 'XLim', [400 600])
    format_plot
    set(gca, 'PlotBoxAspectRatio', [2 1 1])

    subplot(212)
    edges = bw(ii)-80:2:bw(ii)+80;
    NN = histc(bwprime(ii,:), edges);
    bar(edges, NN,'histc')
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','none')
%     set(gca, 'YLim', [0 600])
%     set(gca, 'XLim', [0 200])
    format_plot
    set(gca, 'PlotBoxAspectRatio', [2 1 1])
end

%%
%             data = iddata(win.*curSegment,[],1); % Package input
%             m = armax(data,[lpcOrder zOrder]); % Call estimator with desired model orders
% 
%             allCoeffsP(i,:) = m.a;
%             allCoeffsZ(i,:) = m.c;
