addpath(genpath('../')); % Paths

f = [500];
bw = [100];
Z= []; Zbw = [];
fs = 8000;
numTrials = 10000;

[N, D] = fb2tf(f, bw, Z, Zbw, fs);
 
% figure, freqz(1, D, 100, fs)

%%
fprime = zeros(length(f), numTrials);
bwprime = zeros(length(bw), numTrials);

randn('state',sum(100*clock)); % Seeds

for ii = 1:numTrials
    s = filter(1, D, randn(1, 1000));

    aprime = arcov(s, length(f)*2);

    [z, p, k] = tf2zp(1, aprime);
    kk = find(imag(p) == 0);
    p(kk) = [];
    p = p(1:2:end);

    fp = angle(p)/(2*pi)*fs;
    bwp = -log(abs(p))*fs/pi;
    
    [fp, jj] = sort(fp);
    bwp = bwp(jj);
    
    fprime(:, ii) = fp;
    bwprime(:, ii) = bwp;
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
    disp(['Freq RMSE: ', num2str(sqrt(mean((fprime(ii,:)-f(ii)).^2)))])
    disp(['BW RMSE: ', num2str(sqrt(mean((bwprime(ii,:)-bw(ii)).^2)))])
    
    disp(' ')
end

format_plot
ylim([0 600])

%%
figure

subplot(211)
edges = 400:2:600;
NN = histc(fprime, edges);
bar(edges, NN,'histc')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','none')
set(gca, 'YLim', [0 1200])
set(gca, 'XLim', [400 600])
format_plot
set(gca, 'PlotBoxAspectRatio', [2 1 1])

subplot(212)
edges = 20:2:200;
NN = histc(bwprime, edges);
bar(edges, NN,'histc')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','none')
set(gca, 'YLim', [0 600])
set(gca, 'XLim', [0 200])
format_plot
set(gca, 'PlotBoxAspectRatio', [2 1 1])

%%
%             data = iddata(win.*curSegment,[],1); % Package input
%             m = armax(data,[lpcOrder zOrder]); % Call estimator with desired model orders
% 
%             allCoeffsP(i,:) = m.a;
%             allCoeffsZ(i,:) = m.c;
