function rate = computeRate_RealAwgnChan_MassPoints(varianceRealNoise, symbols, symbolProbabilities)
%computeRate_RealAwgnChan_MassPoints: Computes the achievable information
% rate over a real-valued peak-power-limited AWGN channel, given the parameters
% of the mass points that describe a discrete  transmit distribution.
%
%   Usage: rate = computeRate_RealAwgnChan_MassPoints(varianceRealNoise, symbols, symbolProbabilities)
%
%   Input varianceRealNoise is the AWGN variance.
%   Input symbols is a vector containing the mass point locations.
%   Input symbolProbabilities is a vector containing the mass points
%   probabilities.
%   Output rate is the information rate in bit per channel use (bpcu).
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

   if ~exist('symbolProbabilities','var')
      symbolProbabilities = ones(size(symbols)) / length(symbols);
   end
   
   %if abs(sum(symbolProbabilities)-1) > eps
   %   warning('computeRate_ComplexAwgnChan_MassPoints: Argument symbolProbabilities does not sum to 1')
   %end

   % shorthand aliases:
   sm = symbols;
   pm = symbolProbabilities;
   sigma = sqrt(varianceRealNoise);

   integrand = @(y) computeRate_RealAwgnChan_MassPoints_SpecificIntegrand(y, varianceRealNoise, sm, pm);
   lb = min(sm) - 5*sigma;
   ub = max(sm) + 5*sigma;
   hY = integral(integrand, lb, ub);
   rate = hY - .5 * log2(2*pi*exp(1)*varianceRealNoise);
end

function value = computeRate_RealAwgnChan_MassPoints_SpecificIntegrand(y, varOfRealNoise, sm, pm)
   fY = computeRate_RealAwgnChan_MassPoints_PDFfy(y, varOfRealNoise, sm, pm);
   value = -fY .* log2(fY);
   value(isinf(value)) = 0;
   value(isnan(value)) = 0;
end

function fY = computeRate_RealAwgnChan_MassPoints_PDFfy(y, varOfRealNoise, sm, pm)
   % Compute the Gaussian mixture:
   fY = zeros(size(y));
   for m = 1 : length(sm)
      d = y - sm(m);
      fY = fY + pm(m)/sqrt(2*pi*varOfRealNoise) * exp(-d.^2 / (2*varOfRealNoise));
   end
end