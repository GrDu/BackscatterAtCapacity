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

SNR_dB = 15;
SNR = 10.^(SNR_dB / 10);

% solution with double integral
resolution = 501;
bVal = linspace(0, 1 + 5/sqrt(SNR), resolution);
db = (max(bVal) - min(bVal)) / (resolution-1);
fb_over_b = nan(1,resolution);
for nb = 1 : resolution
   b = bVal(nb);
   kernel_a = @(a) 2*SNR * 2*a .* exp(-SNR*(a - b).^2) .* besseli(0, 2*SNR*(a.*b), 1); % for scaling details see: help besseli
   fb_over_b(nb) = integral(kernel_a, 0, 1);
end
h_y = log2(2*pi) - sum(bVal .* fb_over_b .* log2(fb_over_b)) * db;
R_ME_DoubleIntegral = log2(SNR/(pi*exp(1))) + h_y

fb = bVal .* fb_over_b;

% plot it
figure(7777), clf
plot(bVal, fb, 'b', 'LineWidth', 1), hold on

% try out the analytical solution for f_b from Shamai (33)
fb_analytical = nan(size(fb));
for nb = 1 : resolution
  b = bVal(nb);
  fb_analytical(nb) = 2*b * (1 - marcumq(b*sqrt(2*SNR),sqrt(2*SNR)));
end

% plot it
plot(bVal, fb_analytical, 'k--', 'LineWidth', 2)

% this seems to work! calculate rate directly with that:
fb_over_b = @(b) 2*(1 - marcumq(b*sqrt(2*SNR),sqrt(2*SNR)));
kernel_b = @(b) b .* fb_over_b(b) .* log2(fb_over_b(b));
h_y = log2(2*pi) - integral(kernel_b, 0, 1 + 5/sqrt(SNR));
R_ME = log2(SNR/(pi*exp(1))) + h_y
