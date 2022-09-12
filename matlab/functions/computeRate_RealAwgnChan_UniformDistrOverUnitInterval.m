function rate = computeRate_RealAwgnChan_UniformDistrOverUnitInterval(varianceRealNoise)
%computeRate_RealAwgnChan_UniformDistrOverUnitInterval: Computes the
%achievable rate resulting from a transmit signal with uniform distribution
%over the interval [-1,1], at a specified noise variance.
%
%   Usage: rate = computeRate_RealAwgnChan_UniformDistrOverUnitInterval(varianceRealNoise)
%
%   Input varianceRealNoise has to be a positive real number in linear scale (not in dB).
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

   integrand = @(y) computeRate_RealAwgnChan_UDoUI_SpecificIntegrand(y, varianceRealNoise);
   limit = 1 + 5*sqrt(varianceRealNoise);
   hY = integral(integrand, -limit, limit);
   rate = hY - .5 * log2(2*pi*exp(1)*varianceRealNoise);
end

function value = computeRate_RealAwgnChan_UDoUI_SpecificIntegrand(y, varianceRealNoise)
   fY = computeRate_RealAwgnChan_UDoUI_PDFfy(y, varianceRealNoise);
   value = -fY .* log2(fY);
   value(isinf(value)) = 0;
   value(isnan(value)) = 0;
end

function fY = computeRate_RealAwgnChan_UDoUI_PDFfy(y, varianceRealNoise)
   sigma = sqrt(varianceRealNoise); % standard devation
   fY = .5*(qfunc((y-1)/sigma) - qfunc((y+1)/sigma));
end