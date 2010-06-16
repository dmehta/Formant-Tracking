%% plot frequency vs frame with confidence intervals
% repackage
x_estPerFreq = zeros(numTrials, numObs, numStates);
x_errVarPerFreq = zeros(numTrials, numObs, numStates);

for kk = 1:numStates
    for jj = 1:numTrials
        x_estPerFreq(jj, :, kk) = x_est{jj}(kk, :);
        x_errVarPerFreq(jj, :, kk) = x_errVar{jj}(kk, kk, :);
    end
end

obs = 1:numObs;
means = zeros(numStates, numObs);
vOfm = zeros(numStates, numObs);
mOfv = zeros(numStates, numObs);

figure(42), box off, hold on
figure(43), box off, hold on
for kk = 1:numStates
    figure(42)
    [L U ave v] = findCI(x_estPerFreq(:,:,kk), 95);
    fill([obs obs(end:-1:1)], [L U(end:-1:1)], [0.9 0.9 0.9], 'EdgeColor', 'none')
    plot(obs, ave, 'b-', 'MarkerFace', 'b', 'MarkerSize', 1, 'LineWidth', 1)

    means(kk,:) = ave;
    vOfm(kk,:) = v;
    mOfv(kk,:) = mean(x_errVarPerFreq(:,:,kk), 1);
    
    figure(43)
    plot(vOfm(kk,:))
    plot(mOfv(kk,:), 'r')
end

figure(42)
xlabel('Frame')
ylabel('Frequency (Hz)')

figure(43)
xlabel('Frame')
ylabel('Variance (Hz^2)')
legend('Variance of means', 'Mean of variances')