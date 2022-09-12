% 
%--------------------------------------------------------------------------
%   This code repository accompanies the following academic publication:
%   Gregor Dumphart, Johannes Sager, Armin Wittneben,
%   "Load Modulation for Backscatter Communication: Channel Capacity and
%   Capacity-Approaching Finite Constellations."
%   arXiv preprint arXiv:2207.08100, July 2022.  
%   Available online at (open access): https://arxiv.org/abs/2207.08100
%
%   Coded by: Gregor Dumphart (gregord.research@gmail.com), July 2022,
%   ETH Zurich (D-ITET, Wireless Communications Group), Switzerland.
%   BSD 3-Clause License applies. Copyright (c) 2022, Gregor Dumphart.
%--------------------------------------------------------------------------   
clear all
addpath 'functions'

SNR_dB = 15;
SNR = 10^(SNR_dB/10);
[S, pmf, radii, circleIdxOfSymbols] = getConstellationAPSK_OptimizedGivenSNR(SNR);
M = length(S)
K = length(radii)

figure(7), clf
subplot(1,2,1)

phi = linspace(-pi, pi, 181);
for k = 1 : length(radii)
    s = radii(k) * exp(1i * phi); 
    plot(real(s), imag(s), '-', 'Color', .85*[1 1 1], 'LineWidth', .5, 'HandleVisibility', 'off'), hold on
end

plot(real(S), imag(S), 'k.', 'MarkerSize', 12, 'HandleVisibility', 'off')
grid on
axis equal
xlabel('Re \Gamma', 'interpreter', 'tex')
ylabel('Im \Gamma', 'interpreter', 'tex')
xlim([-1.1 1.1])
ylim([-1.4 1.4])

z = (1 + S) ./ (1 - S);
z(abs(z) > 100) = nan;

subplot(1,2,2)

for k = 1 : K
    s = radii(k) * exp(1i * phi);
    zPts = (1 + s) ./ (1 - s);
    zPts(abs(zPts) > 100) = nan;
    plot(real(zPts), imag(zPts), '-', 'Color', .85*[1 1 1], 'LineWidth', .5, 'HandleVisibility', 'off'), hold on
end

plot(real(z), imag(z), 'k.', 'MarkerSize', 12, 'HandleVisibility', 'off')
grid on
axis equal

xlabel('Re z', 'interpreter', 'tex')
ylabel('Im z', 'interpreter', 'tex')
xlim([-.3 2.3]), ylim([-2.4 2.4])
set(gca, 'XTick',  0 : 2)
set(gca, 'YTick', -1 : 1)
