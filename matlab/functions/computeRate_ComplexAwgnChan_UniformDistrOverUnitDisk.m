function rate = computeRate_ComplexAwgnChan_UniformDistrOverUnitDisk(varianceComplexAWGN)
%computeRate_ComplexAwgnChan_UniformDistrOverUnitDisk: Computes the
%achievable rate resulting from a transmit signal with uniform distribution
%over the unit disk, at a specified SNR.
%
%   Usage: rate = computeRate_ComplexAwgnChan_UniformDistrOverUnitDisk(varianceComplexAWGN)
%
%   Input varianceComplexAWGN is the variance of the complex-valued AWGN.
%   It must be a positive real number in linear scale (not in dB).
%   It should not be below 1E-8, although the function is very robust.
%
%   Output rate is the mutual information in bit per channel use (bpcu).
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

   SNR = 1 / varianceComplexAWGN;

   % Describe the PDF of the received-signal amplitude with a formula from the
   % paper Section III-C, which is based on [ShamaiBarDavidTIT1995, Eq.(33)].
   % It uses the Marcum Q function. Then compute the mutual information with
   % equation (19) from the paper, which requires numerical integration:
   fb_over_b = @(b) 2*(1 - marcumq(b*sqrt(2*SNR),sqrt(2*SNR)));
   integrand = @(b) b .* fb_over_b(b) .* log2(fb_over_b(b));
   rate = log2(2*SNR/exp(1)) - integral(integrand, 0, 1 + 5/sqrt(SNR));
   
   % % [For Reference] Slower alternative solution in terms of a double integral
   % % (this doesn't use the analytical Marcum Q function description):
   % resolution = 501;
   % bVal = linspace(0, 1 + 5/sqrt(SNR), resolution);
   % db = (max(bVal) - min(bVal)) / (resolution-1);
   % fb_over_b = nan(1,resolution);
   % for nb = 1 : resolution
   %    b = bVal(nb);
   %    integrand_a = @(a) 2*SNR * 2*a .* exp(-SNR*(a - b).^2) .* ...
   %       besseli(0, 2*SNR*(a.*b), 1); % for scaling details see: help besseli
   %    fb_over_b(nb) = integral(integrand_a, 0, 1);
   % end
   % rate = log2(2*SNR/exp(1)) - sum(bVal .* fb_over_b .* log2(fb_over_b)) * db;
end