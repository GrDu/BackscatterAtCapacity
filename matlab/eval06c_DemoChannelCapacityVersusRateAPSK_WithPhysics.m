% First, the signal-to-noise ratio (SNR) is calculated from physical
% quantities for a setup of an RF source, a backscatter tag, and a receiver.
% They are all on the same line at certain specified distances. They are
% equipped with halfwave dipole antennas whose main lobes are pointing
% along the aforementioned line.
% Then, the script demonstrates the construction, use, and meaning of the
% achievable information rate with a near-capacity-achieving finite PSK or
% APSK constellation, used over a quadrature peak-power-limited AWGN channel.
% 
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

% Propagation model:
%
%   ,------------------------>\
% Source                       + Receiver
%   `---------------> Tag --->/
%
% Let's assume that the locations of the source, tag, and receiver are on
% the same line, in order to keep things simple.

% Geometry parameters:
distanceRxSrc  = 30; % [m]
distanceRxTag  = 2; % [m]
distanceTagSrc = abs(distanceRxSrc - distanceRxTag);

% System parameters:
PSrc = 10E-3; % [W] transmit power of the RF source
fc = 2.4E9;   % [Hz] carrier frequency
fSymbol = 1E6; % [Hz] symbol rate
B = fSymbol; % [Hz] receive filter bandwidth
T = 500; % [K] antenna temperature

% Don't ever change these:
RLambda2 = 73.1; % [Ohm] resistance of a halfwave dipole antenna. 
                 % See Balanis equation (4-93) or
                 % https://en.wikipedia.org/wiki/Dipole_antenna#Half-wave_dipole
c = 299792458; % [m/s] vacuum propagation speed
kBoltzmann = 1.38E-23; % [J/K] Boltzmann constant

% Antenna resistances
RSrc = RLambda2; % [Ohm] source
RTag = RLambda2; % [Ohm] tag
RRx  = RLambda2; % [Ohm] receiver

% Derived mutual impedances between antennas:
% Let's assume half-wave dipole antennas. The underlying theory is given in
% the book by Balanis, "Antenna Theory: Analysis and Design", in Cpt. 4.6,
% together with terms for the near- and mid-field and some circuit theory
% (e.g., see PhD thesis Gregor).
k = 2*pi*fc/c; % [rad/m] wavenumber
G = 1.643; % maximum directivitiy, i.e. antenna gain, see Balanis equation (4-91)
kr = k * distanceRxSrc; % shorthand notation
ZRxSrc = G * RLambda2 * (1/kr - 1i/kr^2 - 1/kr^3) * exp(-1i*kr); % [Ohm] mutual impedance
kr = k * distanceRxTag;
ZRxTag = G * RLambda2 * (1/kr - 1i/kr^2 - 1/kr^3) * exp(-1i*kr); % [Ohm] mutual impedance
kr = k * distanceTagSrc;
ZTagSrc= G * RLambda2 * (1/kr - 1i/kr^2 - 1/kr^3) * exp(-1i*kr); % [Ohm] mutual impedance

% Derived complex phasors:
iSrc = sqrt(PSrc / RSrc) * exp(2i*pi*rand(1)); 
                          % [A] source antenna feed current with random phase
vTagInd = ZTagSrc * iSrc; % [V] induced voltage at the tag
vRxInd  = ZRxSrc  * iSrc; % [V] direct-path induced voltage

% Derived AWGN quantities:
N0 = kBoltzmann * T; % [J = W/Hz] noise spectral density
sigma_vN = sqrt(4 * N0 * B * RRx); % [V] standard deviation of receiver noise voltage vN
i0 = vTagInd / (2 * RTag); % [A] current describing the tag antenna feed current "disk"
v0 = ZRxTag * i0;
SNR = abs(v0)^2 / sigma_vN^2; % signal-to-noise ratio (rho)
SNR_PerDim = 2*SNR; % SNR per dimension (real- or imaginary part)
SIR = abs(v0)^2 / abs(vRxInd)^2; % signal-to-interference ratio
SNR_dB = 10*log10(SNR) % [dB]
SIR_dB = 10*log10(SIR) % [dB]
sigma_w = 1 / sqrt(SNR);
sigma_w_PerDim = 1 / sqrt(SNR_PerDim);

% Modulation parameters:
N = 15E3; % coding block length (ToDo: set at discretion, define preamble, ...)
doPurelyReactiveLoadModulation = false;

% Create the information-theoretic reference: what's the channel capacity
% at this SNR?
if doPurelyReactiveLoadModulation
   RateChannelCapacity = computeRate_ComplexAwgnChan_ConcentricCircles(1/SNR, 1, 1)
else
   [RateChannelCapacity, KOpt, akOpt, qkOpt] = channelCapacity_ComplexAwgnChan_PeakPowerConstr(1,1/SNR)
   RateChannelCapacity % [bpcu] channel capacity
   KOpt % associated number of circles
end
MMinReq = 2^RateChannelCapacity % minimum required constellation size that prevents that
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

iTag = (1 - Gamma) * i0;  % [A] tag current by equation (4)
vN = (randn(1,N) + 1i * randn(1,N)) * sigma_vN/sqrt(2); % [V] noise voltage sequence
%w = (randn(1,N) + 1i * randn(1,N)) / sqrt(2*SNR);

% NOW: Received voltage sequence by equation (6):
v = -ZRxTag * iTag + vN + vRxInd;

% NOW: Receive Processing

% % To prevent stupid errors, delete all quantities that the receiver can't know:
% % (however, you might regret this when doing evaluations afterwards)
% clear Gamma iTag vTagInd vRxInd v0 i0 iSrc PSrc kr sigma_vN ...
%       ZRxSrc ZRxTag ZTagSrc SNR SNR_PerDim SIR SNR_dB sigma_w sigma_w_PerDim

% Estimate the direct path:
% Here we just use the true value (which an actual receiver would not know).
vRxInd_estimate = vRxInd;

% Compensate the direct path according to equation (7):
vDiff = v - vRxInd_estimate;

% Channel estimation: Estimate a reference voltage phasor that a noiseless and
% interference-free receiver would experience when iT = i0 (i.e. x = 0, z = 1)
% is transmitted.
% Here we just use the true value (which an actual receiver would not know).
v0_estimate = v0;

% Compensate the channel phase shift and apply automatic gain control (AGC)
% according to equation (9):
y = 1 + vDiff / v0_estimate;

% NOW: Evaluation

% Plot the transmitted symbols versus received samples:
figure(1), clf
plot(real(y), imag(y), 'r.'), hold on
plot(real(S), imag(S), 'kx', 'MarkerSize', 9, 'LineWidth', 3)
xlabel('Re(y)'), ylabel('Im(y)'), axis equal
a = 1 + 5*sigma_w_PerDim; xlim([-a a]), ylim([-a a])

% Quantities that might be useful for the channel encoder/decoder:
EuclideanDistances = abs(S.' - S);
PairwiseErrorProbs = qfunc(EuclideanDistances / sigma_w_PerDim);
PairwiseErrorProbs(1 : M+1 : end) = nan; % eliminate diagonal

% When computing LLRs, don't forget that:
% *) The symbol probabilities of the capacity-approaching APSK modulation 
%    are non-uniform!
% *) The variance of the complex-valued AWGN w is sigma_w^2 = 1/SNR. 
% *) The variance per dimension (real- and imaginary part of w) is 1/(2*SNR).
% *) The standard deviation per dimension is sigma_w_PerDim = 1/sqrt(2*SNR).


