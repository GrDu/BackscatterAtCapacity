function [S, radii, circleIdxOfSymbols, pmf] = getConstellationAPSK(M)
%getConstellationAPSK:
% Returns an APSK constellation that fits into the unit disk.
% Canonical choices are made for the circle radii.
% Uniform symbol probabilities are chosen.
% 
%   Inputs:
%   M is the constellation size. It has to be positive integer that is either
%       a square number or smaller than 9 (in which case a PSK is returned).
%
%   Outputs:
%   S is a vector with the complex-valued symbols.
%   radii is a vector with the different radii (i.e. absolute values) of the
%       symbols.
%   circleIdxOfSymbols is a vector containing indices that establish the
%       symbol-to-radius map, which is often useful
%   pmf is a vector with the probability mass function of the symbols.
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

if M <= 8
    warning('getConstellationAPSK: M is small. Returning a PSK instead.')
    S = getConstellationPSK(M);
    radii = [1];
    circleIdxOfSymbols = ones(1,M);
    pmf = ones(1,M) / M;
    return
end

% Construct APSK by transforming a QAM according to paper
% "QAM to Circular Isomorphic Constellations" by Kayhan, IEEE 2016.
S = getConstellationQAM(M);
for m = 1 : M
    s = S(m);
    if s ~= 0
       S(m) = sqrt(2) * max(abs([real(s),imag(s)])) * s / abs(s);
    end
end
S = S / max(abs(S)); % only done for additional numerical safety

% % APSK acoording to DVB-S2X standard,
% % https://ch.mathworks.com/help/comm/ref/dvbsapskmod.html
% S = dvbsapskmod(0 : M-1, M, 's2x');

% Some rudimentary pre-ordering
[~,idx] = sort(abs(S) - 1E-6 * angle(S),'descend');
S = S(idx);

% Figure out the radii of the APSK constellation
[radii,~,circleIdxOfSymbols] = uniquetol(abs(S),1E-6);
K = length(radii);

% Flip the order (for formal consistency with the paper)
radii = flip(radii);
circleIdxOfSymbols = K - circleIdxOfSymbols + 1;

% Fix non-equidistant angular spacing on a circle and introduce a smart
% rotation that improves Euclidean distances
for k = 1 : length(radii)
    idx = (circleIdxOfSymbols == k);
    M_k = sum(idx);
    if M_k == 1
       S(idx) = 0;
    else
       anglesProper = 2*pi * ((0 : M_k-1)' + 1/2 * mod(k,2)) / M_k;
       S(idx) = radii(k) * exp(1i*anglesProper);
    end
end

% Set the symbol pmf to a uniform distribution
pmf = 1/M * ones(size(S));