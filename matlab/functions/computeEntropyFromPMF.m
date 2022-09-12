function H = computeEntropyFromPMF(pmf)
%computeEntropyFromPMF: Computes the entropy H of a discrete random variable.
%
%   Usage: H = computeEntropyFromPMF(pmf).
%
%   Inputs:
%   pmf is the probability mass function of the RV. It must either be a
%       vector or a matrix. If pmf is a vector, it should fulfill
%       sum(pmf) == 1, which is not checked by this function (because
%       secondary uses are possible, cf. computeRate_ComplexAwgnChan_MassPoints).
%       If pmf is a matrix, then the entropy is calculated row-wise and the
%       sum of each row of pmf should be 1.
%
%   Output:
%   H is the calculated entropy. If the input pmf was a vector then H is a
%       scalar. If pmf was a matrix then H is a column vector.
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

% if abs(sum(pmf)-1) > eps
%     warning('computeEntropyFromPMF: Argument pmf does not sum to 1')
% end

tmp = -pmf .* log2(pmf);
tmp(isinf(tmp)) = 0;
tmp(isnan(tmp)) = 0;
if size(tmp,1) > 1 && size(tmp,2) > 1 % oh, it's a matrix!
    H = sum(tmp, 2);
else
    H = sum(tmp);
end