% Demonstrates the construction, use, and meaning of the achievable
% information rate with a finite PSK or APSK constellation, used over a
% quadrature peak-power-limited AWGN channel.
%
% The mathematics in here are guided by the paper stated right below.
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
clear all
addpath 'functions'

% Tune these parameters at will
SNR_dB = 15;
N = 15E3; % coding block length (set at discretion)
doPurelyReactiveLoadModulation = false;

% Derived parameters
SNR = 10^(SNR_dB/10);
SNR_PerDim = 2*SNR; % SNR per dimension (real- or imaginary part)
sigma_w = 1 / sqrt(SNR);
sigma_w_PerDim = 1 / sqrt(SNR_PerDim);

% Create the information-theoretic reference: what's the channel capacity
% at this SNR?
if doPurelyReactiveLoadModulation
   RateChannelCapacity = computeRate_ComplexAwgnChan_ConcentricCircles(1/SNR, 1, 1)
else
   [RateChannelCapacity, KOpt, akOpt, qkOpt] = channelCapacity_ComplexAwgnChan_PeakPowerConstr(1,1/SNR)
   RateChannelCapacity % [bpcu] channel capacity
   KOpt % associated number of circles
end
MMinReq = ceil(2^RateChannelCapacity) % minimum required constellation size that prevents that
                                      % the information rate will definitely be bottlenecked by
                                      % the discretization of the transmit distribution.

% NOW: Creating the  transmit signal
if doPurelyReactiveLoadModulation
   margin = 1; % number of extra uncoded bpcu to reduce discretization loss
   M = 2^(ceil(RateChannelCapacity)+margin)
   S = getConstellationPSK(M);
   symbolProbs = ones(1,M) / M
else
    % % APSK with specified M. Number of circles K is constructed from M.
    % [S, symbolProbs] = getConstellationAPSK_OptimizedGivenM(16);
    
    % % APSK with specified M. Optional parameter SNR: try to optimize
    % % for target SNR (while keeping the specified M):
    % [S, symbolProbs] = getConstellationAPSK_OptimizedGivenM(16,SNR);
    
    % APSK optimized for specified SNR. M and K result from optimization:
    [S, symbolProbs] = getConstellationAPSK_OptimizedGivenSNR(SNR);
    
    M = length(S)
    symbolProbs
end
H_source = computeEntropyFromPMF(symbolProbs) % this is an upper bound for the achievable rate

bpcuUncoded = log2(M)
if bpcuUncoded < RateChannelCapacity
   % This will not be triggered but it's important to have that aspect in mind!
   warning('Information rate will be bottlenecked by the discretization of the transmit distribution!')
end

% More information theory: What is the achievable information rate with
% this specific symbol alphabet (together with the symbol probabilities)?
RateWithAlphabet = computeRate_ComplexAwgnChan_MassPoints(1/SNR, S, symbolProbs)
discretizationLoss_bpcu = RateChannelCapacity - RateWithAlphabet

% Sample the symbols from alphabet S according to the symbol probabilities:
Gamma = randsample(S,N,true,symbolProbs);

% Noisy received signal y:
w = (randn(1,N) + 1i * randn(1,N)) * sigma_w_PerDim;
y = Gamma + w;

% NOW: Evaluation

% Plot the transmitted symbols versus received samples:
figure(60201), clf
plot(real(y), imag(y), 'r.'), hold on
plot(real(S), imag(S), 'kx', 'MarkerSize', 9, 'LineWidth', 3)
xlabel('Re(y)'), ylabel('Im(y)'), axis equal
a = 1 + 5*sigma_w_PerDim; xlim([-a a]), ylim([-a a])

% Quantities that might be useful for the channel encoder/decoder:
EuclideanDistances = abs(S.' - S);
PairwiseErrorProbabilities = qfunc(EuclideanDistances / sigma_w_PerDim);
PairwiseErrorProbabilities(1 : M+1 : end) = nan; % eliminate diagonal

% When computing LLRs, don't forget that:
% *) The symbol probabilities of the capacity-approaching APSK modulation 
%    are non-uniform!
% *) The variance of the complex-valued AWGN w is sigma_w^2 = 1/SNR. 
% *) The variance per dimension (real- and imaginary part of w) is 1/(2*SNR).
% *) The standard deviation per dimension is sigma_w_PerDim = 1/sqrt(2*SNR).
