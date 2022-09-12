% Computes interesting information rates versus the signal-to-noise ration
% (SNR), using the appropriate function calls. Then plots the results.
%
% This generates Fig.3 (with some add-ons) of the paper stated below.
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

SNR_dB_min = -30;
SNR_dB_max = 50;

SNR_dB = SNR_dB_min : 1/3 : SNR_dB_max;
SNR = 10.^(SNR_dB / 10);
varianceComplexNoise = 1 ./ SNR;
varianceRealNoise = varianceComplexNoise / 2;

[R_max, R_UD, R_react, R_resist, R_UD_real] = deal(nan(size(SNR)));

for n = 1 : length(SNR)
   fprintf('Progress: SNR = %.3f dB\n', SNR_dB(n));
   R_max(n) = channelCapacity_ComplexAwgnChan_PeakPowerConstr(1,varianceComplexNoise(n));
   R_UD(n) = computeRate_ComplexAwgnChan_UniformDistrOverUnitDisk(varianceComplexNoise(n));
   R_react(n) = computeRate_ComplexAwgnChan_ConcentricCircles(varianceComplexNoise(n), 1, 1);
   R_resist(n) = channelCapacity_RealAwgnChan_PeakPowerConstr(1,varianceRealNoise(n));
   R_UD_real(n) = computeRate_RealAwgnChan_UniformDistrOverUnitInterval(varianceRealNoise(n));
end
R_UBX = SNR * log2(exp(1));
R_UB  = log2(1 + SNR);
R_LB  = log2(1 + SNR/exp(1));
R_asymp_react  = .5 * log2(4*pi*SNR/exp(1));
R_asymp_resist = .5 * log2(1 + 4*SNR/(pi*exp(1)));

figure(2), clf

subplot(1,2,1)

plot(SNR_dB, R_UBX, 'r--', 'LineWidth', 2, 'DisplayName', 'upper bound \rho \cdot log_2(e)'), hold on
plot(SNR_dB, R_UB , 'r:' , 'LineWidth', 2, 'DisplayName', 'upper bound log_2(1 + \rho)'), hold on
plot(SNR_dB, R_max, 'b-' , 'LineWidth', 2, 'DisplayName', 'channel capacity, general passive load'), hold on
plot(SNR_dB, R_UD , 'b-.', 'LineWidth', 2, 'DisplayName', 'ach. rate, unif. distr. on complex unit disk')
plot(SNR_dB, R_LB , 'b:' , 'LineWidth', 2, 'DisplayName', 'lower bound log_2(1 + \rho/e)')
plot(SNR_dB, R_react,  'k--', 'LineWidth', 2, 'DisplayName', 'channel capacity, purely reactive load')
plot(SNR_dB, R_asymp_react, 'k:', 'LineWidth', 2, 'DisplayName', 'reactive asymptote .5 log_2(4\pi\rho/e)')
plot(SNR_dB, R_resist,  'g--', 'LineWidth', 2, 'DisplayName', 'channel capacity, purely resistive load')
plot(SNR_dB, R_UD_real, 'g-.', 'LineWidth', 2, 'DisplayName', 'ach. rate, unif. distr. on real-valued [-1,1]')
plot(SNR_dB, R_asymp_resist, 'g:', 'LineWidth', 2, 'DisplayName', 'resistive lower bound .5 log_2(1 + 4\rho/(e\pi))')

grid on
xlabel('signal-to-noise ratio \rho [dB]')
ylabel('achievable information rate [bpcu]')
legend('show', 'Location', 'NorthWest')
xlim([min(SNR_dB), max(SNR_dB)]),
ylim([0, ceil(max(R_max))]),
xticks(-100 : 10 : 100)

subplot(1,2,2)

semilogy(SNR_dB, R_UBX, 'r--', 'LineWidth', 2, 'DisplayName', 'upper bound \rho \cdot log_2(e)'), hold on
semilogy(SNR_dB, R_UB , 'r:' , 'LineWidth', 2, 'DisplayName', 'upper bound log_2(1 + \rho)'), hold on
semilogy(SNR_dB, R_max, 'b-' , 'LineWidth', 2, 'DisplayName', 'channel capacity, general passive load'), hold on
semilogy(SNR_dB, R_UD , 'b-.', 'LineWidth', 2, 'DisplayName', 'achievable rate, maximum-entropy signaling')
semilogy(SNR_dB, R_LB , 'b:' , 'LineWidth', 2, 'DisplayName', 'lower bound log_2(1 + \rho/e)')
semilogy(SNR_dB, R_react,  'k--', 'LineWidth', 2, 'DisplayName', 'channel capacity, purely reactive load')
semilogy(SNR_dB, R_asymp_react, 'k:', 'LineWidth', 2, 'DisplayName', 'reactive asymptote .5 log_2(4\pi\rho/e)')
semilogy(SNR_dB, R_resist,  'g--', 'LineWidth', 2, 'DisplayName', 'channel capacity, purely resistive load')
semilogy(SNR_dB, R_UD_real, 'g-.', 'LineWidth', 2, 'DisplayName', 'ach. rate, unif. distr. on real-valued [-1,1]')
semilogy(SNR_dB, R_asymp_resist, 'g:', 'LineWidth', 2, 'DisplayName', 'resistive lower bound .5 log_2(1 + 4\rho/(e\pi))')

grid on
xlabel('signal-to-noise ratio \rho [dB]')
ylabel('achievable information rate [bpcu]')
legend('show', 'Location', 'NorthWest')
xlim([min(SNR_dB), 10]),
ylim([1E-3, 1E2]),
xticks(-100 : 10 : 100)
