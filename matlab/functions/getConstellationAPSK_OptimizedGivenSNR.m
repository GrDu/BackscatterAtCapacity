function [S, pmf, radii, circleIdxOfSymbols] = getConstellationAPSK_OptimizedGivenSNR(DesignSNR)
% For a specified target design signal-to-noise ratio (SNR), an
% amplitude-and-phase-shift-keying (APSK) constellation on the unit disk is
% returned.
%
%   Inputs:
%   DesignSNR is the target design SNR in linear scale (not in dB).
%
%   Outputs:
%   S is a vector with the complex-valued symbols.
%   radii is a vector with the different radii (i.e. absolute values) of the
%       symbols.
%   circleIdxOfSymbols is a vector containing indices that establish the
%       symbol-to-radius map, which is often useful.
%   pmf is a vector with the probability mass function of the symbols.
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

[RateMax, K, ak, qk] = channelCapacity_ComplexAwgnChan_PeakPowerConstr(1,1/DesignSNR);

if isnan(K)
    load 'channelCapacity_ComplexAwgnChan_PeakPowerConstr_4to40dB.mat'
    warning('getConstellationAPSK_OptimizedGivenSNR: DesignSNR = %.2f dB is extremely high. Can not assign optimal distribution params. Using params for %.2f dB instead.', 10*log10(DesignSNR), 10*log10(SNR(end)))
    DesignSNR = SNR(end); % from load data
    [RateMax, K, ak, qk] = channelCapacity_ComplexAwgnChan_PeakPowerConstr(1,1/DesignSNR);
end

% Choose constellation size M (must be a square number) and calculate K:
if length(ak) == 1
    M = 2^ceil(RateMax);
elseif ak(end) > 0 % it's a proper circle -> use even-numbered square number M
    M = 4*K^2; % = sum((4 : 8 : 4+8*(K-1)))
else % ok, it's a point -> use an odd-numbered square number M
    M = 1 + 4*K*(K-1); % = 1 + sum(8 : 8 : 8*(K-1))
    %M = (K+1)^2; 
end

if log2(M) < RateMax % This should never happen :-)
    warning('getConstellationAPSK_OptimizedGivenSNR: Information rate will definitely be bottlenecked by the discretization of the transmit distribution!')
end

% Generate a basic APSK constellation:
[S, radiiBasic, circleIdxOfSymbols, pmf] = getConstellationAPSK(M);
if length(radiiBasic) ~= K % This shouldn't happen either, by design :-)
    error('getConstellationAPSK_OptimizedGivenSNR: Mismatching number of radii!')
end

% Set radii, angles, and constellation PMF appropriately
for k = 1 : length(radiiBasic)
    idx = (circleIdxOfSymbols == k);
    if radiiBasic(k) > 0
        S(idx) = S(idx) * ak(k) / radiiBasic(k); 
    end
    pmf(idx) = qk(k) / sum(idx);
end
radii = ak;