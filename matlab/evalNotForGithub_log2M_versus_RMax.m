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

SNR_dB_min = -10;
SNR_dB_max = 30;

SNR_dB = SNR_dB_min : .05 : SNR_dB_max;
SNR = 10.^(SNR_dB / 10);
varianceComplexNoise = 1 ./ SNR;
varianceRealNoise = varianceComplexNoise / 2;

[R_max, K_evo, M_evo, R_react] = deal(nan(size(SNR)));

for n = 1 : length(SNR)
   [rate, K, ak, qk] = channelCapacity_ComplexAwgnChan_PeakPowerConstr(1,varianceComplexNoise(n));
   R_max(n) = rate;
   K_evo(n) = K;

   if length(ak) == 1
      M_evo(n) = 2^ceil(rate);
   elseif ak(end) > 0 % it's a proper circle -> use even-numbered square number M
      M_evo(n) = 4*K^2; % = sum((4 : 8 : 4+8*(K-1)))
   else % ok, it's a point -> use an odd-numbered square number M
      M_evo(n) = (K+1)^2; % = 1 + sum(8 : 8 : 8*(K-1))
   end

   R_react(n) = computeRate_ComplexAwgnChan_ConcentricCircles(varianceComplexNoise(n), 1, 1);
end

figure(100001), clf

semilogy(SNR_dB, R_max, 'b-' , 'LineWidth', 2, 'DisplayName', 'channel capacity, general passive load'), hold on
semilogy(SNR_dB, log2(M_evo),  'r--', 'LineWidth', 2, 'DisplayName', 'upper bound with chosen discretization')
semilogy(SNR_dB, R_react,  'k--', 'LineWidth', 2, 'DisplayName', 'channel capacity, purely reactive load')

grid on
xlabel('signal-to-noise ratio \rho [dB]')
ylabel('achievable information rate [bpcu]')
legend('show', 'Location', 'NorthWest')
xlim([min(SNR_dB), max(SNR_dB)]),
%ylim([0, ceil(max(R_max))]),
xticks(-100 : 10 : 100)

figure(100002), clf

semilogy(SNR_dB, M_evo)
