function rate = computeRate_ComplexAwgnChan_MassPoints(noiseVariance, symbols, symbolProbabilities, resolution)
%computeRate_ComplexAwgnChan_MassPoints: Computes the achievable information
% rate over a complex-valued peak-power-limited AWGN channel, given the
% parameters of the mass points that describe a discrete transmit distribution.
%
%   Usage:
%   rate = computeRate_ComplexAwgnChan_MassPoints(noiseVariance, symbols, symbolProbabilities, resolution)
%
%   Inputs:
%   noiseVariance is the variance of the complex-valued additive white
%       Gaussian noise (AWGN).
%   symbols is a vector holding the complex-valued constellation
%       points s_m.
%   symbolProbabilities holds the smybol probabilites p_m. This is an
%       optional parameter; a uniform distribution is assumed if omitted.
%   resolution defines the integration resolution per dimension for the 
%       underlying differential entropy computations. This is an optional
%       parameter; a resolution of 120 is assumed if omitted.
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

if ~exist('symbolProbabilities','var')
   symbolProbabilities = ones(size(symbols)) / length(symbols);
end

%if abs(sum(symbolProbabilities)-1) > eps
%   warning('computeRate_ComplexAwgnChan_MassPoints: Argument pmf does not sum to 1')
%end

M = length(symbols);
sigma = sqrt(noiseVariance);
yDiskRadius = max(abs(symbols)) + 4*sigma/sqrt(2); % consider 5 std.dev. (per dimension) excess over unit disk

if ~exist('resolution','var')
   resolution = 120;
end
val = linspace(-yDiskRadius, yDiskRadius, resolution+1);
y = val + 1i * val'; % generate a nice grid
y = y(:); % collapse it into a vector
y = y(abs(y) <= yDiskRadius); % ignore points outside a certain disk
differential = 2*yDiskRadius / resolution;

% Compute the received-signal PDF, a two-dimensional Gaussian mixture:
pdf_Y = zeros(size(y));
for m = 1 : M
   % https://en.wikipedia.org/wiki/Complex_normal_distribution#Probability_density_function
   d = abs(y - symbols(m));
   tmp = 1/(pi*sigma^2) * exp(-d.^2 / sigma^2);
   pdf_Y_given_sm = tmp / sum(tmp) / differential^2; % normalize
   pdf_Y = pdf_Y + pdf_Y_given_sm * symbolProbabilities(m);
end

h_Y = computeEntropyFromPMF(pdf_Y) * differential^2; % differential entropy of received signal
h_W = log2(pi*exp(1)*sigma^2); % differential entropy of the noise
rate = h_Y - h_W; % mutual information; a difference of differential entropies