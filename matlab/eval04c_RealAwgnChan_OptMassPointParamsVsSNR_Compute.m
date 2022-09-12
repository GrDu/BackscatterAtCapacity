% Computes the parameters of the (finite number of) mass points that
% describe the capacity-achieving transmit distribution over a
% peak-power-constrained real-valued AWGN channel. This is done for many
% different signal-to-noise ratio (SNR) values with a thorough sweep.
%
% The underlying theory is given in the paper:
% Joel G. Smith 1971, "The Information Capacity of Amplitude- and
% Variance-Constrained Scalar Gaussian Channels".
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

% % See the paper: Joel G. Smith 1971, "The Information Capacity of
% % Amplitude- and Variance-Constrained Scalar Gaussian Channels".
% A_min = 1;
% A_max = 40;
% SNR_min = .5 * A_min^2;
% SNR_max = .5 * A_max^2;
% %SNR_dB_min = floor(10*log10(SNR_min));
% %SNR_dB_max = ceil(10*log10(SNR_max));

% Config in "the other direction"
SNR_DefPaper2022_dB_min = 0;
SNR_DefPaper2022_dB_max = 40;
%SNR_PaperDef_min = 10^(SNR_PaperDef_dB_min/10);
%SNR_PaperDef_max = 10^(SNR_PaperDef_dB_max/10);
%A_min = sqrt(2 * SNR_PaperDef_min)
%A_max = sqrt(2 * SNR_PaperDef_max)

SNR_DefPaper2022_dB = SNR_DefPaper2022_dB_min : .1 : SNR_DefPaper2022_dB_max;
SNR_DefPaper2022 = 10.^(SNR_DefPaper2022_dB / 10);
varianceRealNoise = 1 ./ (2*SNR_DefPaper2022);
%A = sqrt(2 * SNR); % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

R_max_real = nan(1, length(varianceRealNoise));
MOnesided_evolution = nan(1, length(varianceRealNoise)); % Number of symbols per side
smOnesided_evolution = nan(99, length(varianceRealNoise)); % mass point (symbol) values
pmOnesided_evolution = nan(99, length(varianceRealNoise)); % mass point probabilities
% 99 is just some sufficiently large max. number of mass points per side

% Explanation: "per side" and "onesided" mean that we only save the
% positive-valued symbols because of the symmetry of the optimal transmit
% distribution (and every candidate for optimality). In detail, if sm is a
% mass point, then -sm is a mass point too. The extremal values +1 and -1
% are always mass points. We don't mind a possible duplicate +0 and -0,
% which does not cause any problems (because of how we split probability).

tic
for n = 1 : length(varianceRealNoise)
   if n == 1
      MOnesided = 1;
      smOnesided_init = 1;
      pmOnesided_init = 1;
   else
      MOnesided = MOnesided_evolution(n-1);
      smOnesided_init = smOnesided_evolution(1 : MOnesided, n-1)';
      pmOnesided_init = pmOnesided_evolution(1 : MOnesided, n-1)';
   end
   
   [rate, smOnesided, pmOnesided] = computeMaxRate_RealAwgnChan_MassPoints(varianceRealNoise(n), smOnesided_init, pmOnesided_init);
   
   % now test whether adding another point improves the data rate:
   pMNewOnesided_init = 1E-2 * pmOnesided(end);
   pmOnesided_extra = [(1-pMNewOnesided_init)*pmOnesided, pMNewOnesided_init];
   smOnesided_extra = [smOnesided, 0];
   rate_extraPt = computeMaxRate_RealAwgnChan_MassPoints(varianceRealNoise(n), smOnesided_extra, pmOnesided_extra);
   rateGainFromExtraPt = rate_extraPt - rate;
   
   % Did the added mass point increase the rate appreciably?
   if rateGainFromExtraPt > 1E-6 * rate
      % Okay, keep the extra trial mass point:
      MOnesided = MOnesided + 1;
      rate = rate_extraPt;
      smOnesided = smOnesided_extra;
      pmOnesided = pmOnesided_extra;
   end

   % Reference information rates (at high SNR, these should be just slightly smaller):
   R_UD_real(n) = computeRate_RealAwgnChan_UniformDistrOverUnitInterval(varianceRealNoise(n));
   R_LB_real(n) = .5 * log2(1 + 4*SNR_DefPaper2022(n)/(pi*exp(1)));
   
   % debug output
   fprintf('SNR = %.2f dB, M_ = %d, R_max = %.3f bpcu, R_UD/R_max = %.6f, R_LB/R_max = %.6f, rateGainFromExtraPt = %.6f bpcu\n', ...
      SNR_DefPaper2022_dB(n), MOnesided, rate, R_UD_real(n) / rate, R_LB_real(n) / rate, rateGainFromExtraPt)
   %smOnesided % more debug output
   %pmOnesided % more debug output
   
   % save the result
   MOnesided_evolution(n) = MOnesided;
   smOnesided_evolution(1 : MOnesided, n) = smOnesided';
   pmOnesided_evolution(1 : MOnesided, n) = pmOnesided';
   R_max_real(n) = rate;
end
toc

% Drop the extra-allocated memory
pmOnesided_evolution = pmOnesided_evolution(1 : max(MOnesided_evolution), :);
smOnesided_evolution = smOnesided_evolution(1 : max(MOnesided_evolution), :);

% % Save the relevant computed data to a file:
%save('channelCapacity_RealAwgnChan_PeakPowerConstr_INFO.mat','SNR_DefPaper2022','SNR_DefPaper2022_dB','varianceRealNoise',...
%     'R_max_real','R_UD_real','R_LB_real',...
%     'pmOnesided_evolution','smOnesided_evolution','MOnesided_evolution')
