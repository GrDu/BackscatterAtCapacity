% Generates the intensity plots of the probability density functions
% (in terms of real and imaginary parts) of the following complex-valued
% random variables:
% - The capacity-achieving transmit signal Gamma for a quadrature
%   peak-power-limited AWGN channel at a specified SNR.  In backscatter
%   communication, Gamma is the reflection coefficient of the passive tag
%   load.
% - The normalized impedance z = (1+Gamma)/(1-Gamma) associated with a
%   random reflection coefficient Gamma.
%
% This generates Fig.5 and Fig.12 of the paper stated below.
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
load 'channelCapacity_ComplexAwgnChan_PeakPowerConstr_4to40dB'

SNR_dB_choice = 18; % < 4.8, 8, 12, 15, (18), 24
[~,nChoice] = min(abs(SNR_dB - SNR_dB_choice));
a = ak_evolution(:,nChoice);
a = a(~isnan(a))' / a(1);
NRadii = length(a);
pmf = qk_evolution(1 : NRadii, nChoice)';

canvas = 1.2;
resolution = 60;
phi_s  = linspace(0, 2*pi, resolution); % for circles

Gamma = a' * exp(1i*phi_s);
zL = (1 + Gamma) ./ (1 - Gamma);

zCirclesPDF = nan(NRadii, resolution-1);

xL = linspace(-canvas, canvas, resolution);
xL = xL(1 : end-1);
zCirclesPDF(1,:) = 1/pi * 1./(1 + xL.^2);


circumference = 4*pi*a ./ (1 - a.^2);        
for idxCirc = 2 : NRadii
    distances = abs(zL(idxCirc,1:end-1) - zL(idxCirc,2:end));
    zCirclesPDF(idxCirc,:) = 1/(resolution-1) * 1 ./ distances;
end
tmp = diag(pmf) * zCirclesPDF;
%color_z = 1 - tmp / (pmf(1) / pi); % make max at (0,0) black always
color_z = 1 - tmp / (1 / pi); color_z(color_z < 0) = 0; % consistent across plots for different SNR

figure(5), clf

for idxCirc = 1 : NRadii
    circProb_lengthAdjusted = pmf ./ a;
    %color_Gamma = 1 - circProb_lengthAdjusted / max(circProb_lengthAdjusted); % make outer ring black always
    color_Gamma = 1 - circProb_lengthAdjusted; color_Gamma(color_Gamma < 0) = 0; % consistent across plots for different SNR
    plot(real(Gamma(idxCirc,:)), imag(Gamma(idxCirc,:)), 'Color', color_Gamma(idxCirc)*[1 1 1], 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off'), hold on
end

figure(12), clf

for idxCirc = 1 : NRadii
    if idxCirc == 1
        if a(1) ~= 1, error('Unexpected circle sorting'), end
        
        for n = 1 : length(xL)-1
            plot([0 0], [xL(n) xL(n+1)], 'Color', color_z(idxCirc,n)*[1 1 1], 'LineWidth', 2, 'HandleVisibility', 'off'), hold on
        end
    else
        
        for n = 1 : length(zL)-1
            data = [zL(idxCirc,n), zL(idxCirc,n+1)];
            if all(real(data) > 2*canvas) || all(imag(data) > canvas)  || all(imag(data) < -canvas)
                continue;
            end
            plot(real(data), imag(data), 'Color', color_z(idxCirc,n)*[1 1 1], 'LineStyle', '-', 'LineWidth', 2, 'HandleVisibility', 'off'), hold on
            %drawnow
        end
    end
end

figure(5)

grid on, axis equal
%title('disk boundary circle and inner circles')
xlabel('Re \Gamma')
ylabel('Im \Gamma')
xlim([-1.1 1.1])
ylim([-1.1 1.1])

colormap(flipud(gray))
h = colorbar;
set(h, 'Ticks', [0 1])
set(h, 'TickLabels', {'0','high'})
h = ylabel(h, 'probability density');
set(h, 'Position', [1.3625 0.5000 0])
set(gcf, 'Position', [500 500 384 236])

figure(12)

grid on,  axis equal
%title('correspondinng normalized impedance')
xlabel('normalized load resistance r_L')
ylabel('normalized load reactance x_L')
%xlim([0 2*canvas]-.5)
xlim([0 2*canvas]-.2)
ylim([-canvas canvas])
set(gca, 'XTick', -10 : 10) % integers only
set(gca, 'YTick', -10 : 10) % integers only

colormap(flipud(gray))
h = colorbar;
set(h, 'Ticks', [0 1])
set(h, 'TickLabels', {'0','high'})
h = ylabel(h, 'probability density');
set(h, 'Position', [1.3625 0.5000 0])
set(gcf, 'Position', [885 500 384 236])



% figure(44444), clf
% plot(SNR_dB, K_evolution, 'k'), hold on
% 
% a=1/sqrt(3);
% plot(SNR_dB, sqrt(SNR / 3), 'b--'),
% plot(SNR_dB, ceil(sqrt(SNR / 3)), 'b')

