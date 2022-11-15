function [S, pmf, radii, circleIdxOfSymbols, ChosenDesignSNR] = getConstellationAPSK_OptimizedGivenM(M, TargetSNR)
% For a specified constellation size M, an amplitude-and-phase-shift-keying
% (APSK) constellation on the unit disk is returned. Optionally, one can
% also specify a target signal-to-noise ratio (SNR), in which case an
% effort is made to optimize the size-M constellation for that SNR in terms
% of symbol radii and probabilities.
%
%   Usage:
%   [S, pmf, radii, circleIdxOfSymbols, ChosenDesignSNR] = getConstellationAPSK_OptimizedGivenM(M)
%   or
%   [S, pmf, radii, circleIdxOfSymbols, ChosenDesignSNR] = getConstellationAPSK_OptimizedGivenM(M, TargetSNR)
%
%   Inputs:
%   M is the constellation size. It has to be positive integer that is either
%       a square number or smaller than 9 (in which case a PSK is returned).
%   TargetSNR is the target design signal-to-noise ratio (SNR) in linear
%       scale (not in dB). This parameter is optional.
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

% First, generate a basic APSK constellation:
[S, radii, circleIdxOfSymbols, pmf] = getConstellationAPSK(M);
K = length(radii);

% Load precomputed circle data
load 'channelCapacity_ComplexAwgnChan_PeakPowerConstr_4to40dB.mat'
goodIndices = find(K_evolution == K);
if isempty(goodIndices)
   error('getConstellationAPSK: No data for that amout of circles!')
end

if exist('TargetSNR','var') && SNR(goodIndices(1)) <= TargetSNR && TargetSNR <= SNR(goodIndices(end))
    ChosenDesignSNR = TargetSNR;
    [~, ~, ak, qk] = channelCapacity_ComplexAwgnChan_PeakPowerConstr(1,1/TargetSNR);
else
    if exist('TargetSNR','var') && TargetSNR < SNR(goodIndices(1))
        warning('Number of circles K (resulting from the specified M with the employed construction scheme) is larger than optimal K for the specified TargetSNR')
        [~,tmpIdx] = min(abs(SNR(goodIndices) - TargetSNR));
        idxTargetSNR = goodIndices(tmpIdx);
    elseif exist('TargetSNR','var') && SNR(goodIndices(end)) < TargetSNR
        warning('Number of circles K (resulting from the specified M with the employed construction scheme) is smaller than optimal K for the specified TargetSNR')
        [~,tmpIdx] = min(abs(SNR(goodIndices) - TargetSNR));
        idxTargetSNR = goodIndices(tmpIdx);
    else % No TargetSNR was specified, so we choose one ...
    
        % % Canonical choice: Center of the eligible SNR interval
        % idxTargetSNR = round(mean([goodIndices(1), goodIndices(end)])); 
       
        % Choose the upper end of the eligible SNR interval
        % (this prevents that the APSK constellation has an ugly dense point cluster near zero)
        idxTargetSNR = goodIndices(end);
    end

    % Assign the parameter values
    ChosenDesignSNR = SNR(idxTargetSNR);
    ak = ak_evolution(1:K,idxTargetSNR) / ak_evolution(1,idxTargetSNR);
    qk = qk_evolution(1:K,idxTargetSNR);
end

% Set the radii, angles, and constellation PMF appropriately
for k = 1 : K
   idx = (circleIdxOfSymbols == k);
   if radii(k) > 0
       S(idx) = S(idx) * ak(k) / radii(k); 
   end
   pmf(idx) = qk(k) / sum(idx);
end
