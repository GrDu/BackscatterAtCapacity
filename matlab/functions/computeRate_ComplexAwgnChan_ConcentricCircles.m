function rate = computeRate_ComplexAwgnChan_ConcentricCircles(varianceComplexNoise, ak, qk)
%computeRate_ComplexAwgnChan_ConcentricCircles: returns the achievable
%information rate in bit per channel use (bpcu) over a complex-valued AWGN
%channel, whereby the transmit distribution is given in terms of a set of
%concentric circles with certain radii and probabilities.
%
%   Usage:
%   rate = computeRate_ComplexAwgnChan_ConcentricCircles(varianceComplexNoise, ak, qk)
%
%   Inputs:
%   ak ..................... vector with radius values
%   qk ..................... vector with radius probabilities
%   varianceComplexNoise ... variance of complex-valued AWGN
%
%   Outputs:
%   rate ... mutual information in bit per channel use (bpcu)
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

   SNR = 1 / varianceComplexNoise;
   integrand = @(b) computeRateGivenCircles_SpecificIntegrand(b, varianceComplexNoise, ak, qk);
   limit = 1 + 5/sqrt(SNR);
   rate = log2(2*SNR/exp(1)) - integral(integrand, 0, limit, 'ArrayValued', true); % mutual information
end

function value = computeRateGivenCircles_SpecificIntegrand(b, varianceComplexNoise, ak, qk)
   SNR = 1 / varianceComplexNoise;
   fb_over_b = 2*SNR*sum(qk .* exp(-SNR*(b - ak).^2) .* besseli(0, 2*SNR*b*ak, 1)); % scale = 1, see: help besseli
   value = b * fb_over_b * log2(fb_over_b);
   if isnan(value)
      value = 0;
   end
end