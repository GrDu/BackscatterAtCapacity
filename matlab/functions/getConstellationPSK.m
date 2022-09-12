function S = getConstellationPSK(M)
% Returns a phase-shift-keying (PSK) constellation on the unit circle.
%
%   Usage: S = getConstellationPSK(M)
%
%   Input:
%   M is the desired constellation order.
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

if M == 2
    S = [-1,1];
else
    angles = 2*pi * ((0 : M-1) + 1/2) / M;
    S = exp(1i*angles); 
end