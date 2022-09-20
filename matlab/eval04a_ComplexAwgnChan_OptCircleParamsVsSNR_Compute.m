% Computes the parameters of the (finite number of) concentric circles that
% describe the capacity-achieving transmit distribution over a
% peak-power-constrained complex-valued AWGN channel. This is done for many
% different signal-to-noise ration (SNR) values with a thorough sweep.
% Be prepared for some runtime.
%
% This prepares Fig.4 (with add-ons) of the paper stated below. The results
% are plotted in eval04b_ComplexAwgnChan_OptCircleParamsVsSNR_Plot.m
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

SNR_dB_min = 4;
SNR_dB_max = 40;

SNR_dB = SNR_dB_min : .1 : SNR_dB_max;
SNR = 10.^(SNR_dB / 10);
varianceComplexNoise = 1 ./ SNR;


K_evolution = nan(1, length(SNR));
ak_evolution = nan(99, length(SNR)); % 99 is just some sufficiently large max. number of radii
qk_evolution = nan(99, length(SNR));
[R_max, R_UD, R_LB] = deal(nan(size(SNR)));

rateGainFromExtraCircle = NaN;
tic
for n = 1 : length(SNR)
    if n == 1
        K = 1;
        ak_init = 1; % circle radius
        qk_init = 1; % circle probability
    else
        K = K_evolution(n-1);
        ak_init = ak_evolution(1 : K, n-1)'; % circle radii
        qk_init = qk_evolution(1 : K, n-1)'; % circle probabilities
    end

    [rate, ak, qk] = computeMaxRate_ComplexAwgnChan_ConcentricCircles(varianceComplexNoise(n), ak_init, qk_init);
    
    if SNR(n) >= 3.011 % Below this threshold, the theory guarantees that K=1 is optimal!
       
       % Now test whether adding another circle improves the data rate:
       qKNew_init = 1E-2 * qk(end); % initialize the new circle with a meaningful positive probability
       qk_extraCircle = [(1-qKNew_init)*qk, qKNew_init];
       ak_extraCircle = [ak, 0]; % the added circle is initalized with radius zero
       rate_extraCircle = computeMaxRate_ComplexAwgnChan_ConcentricCircles(varianceComplexNoise(n), ak_extraCircle, qk_extraCircle);
       rateGainFromExtraCircle = rate_extraCircle - rate;

       % Did the added circle increase the rate appreciably?
       if rateGainFromExtraCircle > 1E-6 * rate
           % Okay, keep the extra trial circle:
           K = K + 1;
           rate = rate_extraCircle;
           ak = ak_extraCircle;
           qk = qk_extraCircle;
       end
    end
    
    % Reference information rates (at high SNR, these should be just slightly smaller):
    R_UD(n) = computeRate_ComplexAwgnChan_UniformDistrOverUnitDisk(varianceComplexNoise(n));
    R_LB(n) = log2(1 + SNR(n)/exp(1));

    % debug output
    fprintf('SNR = %.2f dB, K = %d, R_max = %.3f bpcu, R_UD/R_max = %.6f, LB/R_max = %.6f, rateGainFromExtraCircle = %.6f bpcu\n', ...
       SNR_dB(n), K, rate, R_UD(n) / rate, R_LB(n) / rate, rateGainFromExtraCircle)
    %ak                   % more debug output
    %qk                   % more debug output
    %lk = qk ./ (2*pi*ak) % more debug output

    % Save the results:
    K_evolution(n) = K;
    ak_evolution(1 : K, n) = ak';
    qk_evolution(1 : K, n) = qk';
    R_max(n) = rate;
end
toc

% Drop the extra-allocated memory
ak_evolution = ak_evolution(1 : max(K_evolution), :);
qk_evolution = qk_evolution(1 : max(K_evolution), :);

% % Save the relevant computed data to a file:
%save('channelCapacity_ComplexAwgnChan_PeakPowerConstr_INFO.mat','SNR','SNR_dB','varianceComplexNoise',...
%     'R_max','R_UD','R_LB','ak_evolution','qk_evolution','K_evolution')
