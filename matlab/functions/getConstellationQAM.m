function S = getConstellationQAM(M)
% Returns a quadrature-amplitude-modulation (QAM) constellation that fits
% into the complex unit disk.
%
%   Usage: S = getConstellationQAM(M)
%
%   Input:
%   M is the desired constellation order. It must be a square number.
%
%   Output:
%   S is a vector with the complex-valued symbols.
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

R = sqrt(M);
if abs(R - round(R)) > eps % is R not an integer?
    error('getConstellationQAM: input M must be a square number')
end

%S = qammod(0 : M-1, M); % Would require communications toolbox
x = linspace(-1/sqrt(2), 1/sqrt(2), R);
S = 1i * x + x';
S = S(:).';
S = S / max(abs(S)); % for numerical safety