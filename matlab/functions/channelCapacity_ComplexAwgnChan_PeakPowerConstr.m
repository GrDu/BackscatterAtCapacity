function [rate, K, ak, qk] = channelCapacity_ComplexAwgnChan_PeakPowerConstr(peakPower, varianceComplexAWGN)
%channelCapacity_ComplexAwgnChan_PeakPowerConstr: Returns the channel
% capacity of a complex-valued (i.e. quadrature) AWGN channel with
% constrained signal peak-power (but no constraint on the signal
% average-power). It also returns the parameters of the (finite number of)
% concentric circles describing the capacity-achieving transmit distribution.
%
%   Usage: Obtain just the rate with
%   rate = channelCapacity_ComplexAwgnChan_PeakPowerConstr(SNR)
%   or also the rate and circle parameters with
%   [rate, K, ak, qk] = channelCapacity_ComplexAwgnChan_PeakPowerConstr(SNR)
%   
%   Inputs:
%   peakPower is the maximum peak power in normal scale (not in dB).
%   varianceComplexAWGN is the variance of the complex-valued
%       additive white Gaussian noise (AWGN). Please note that some sources
%       instead consider noise variance per dimension (real or imaginary
%       part), i.e. half the value.
%   rate is the channel capacity in bit per channel use (bpcu), which is
%       referred to as R_max in the paper.
%   K is the number of circles that describe the capacity-achieving
%       transmit distribution.
%   ak is a length-K vector holding the circle radii.
%   qk is a length-K vector holding the circle probabilities.
%
%   This function relies on loading and interpolating precomputed data.
%   Computing everything from scratch would be quite an ordeal because it
%   requires a fine-grained iteration from low SNR to the specified SNR.
%   This process is implemented in the script
%   eval04a_ComplexAwgnChan_OptCircleParamsVsSNR_Compute.m
%   which is able to generate the employed data.
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

if ~isscalar(peakPower) || ~isreal(peakPower) || peakPower < 0
    error('channelCapacity_ComplexAwgnChan_PeakPowerConstr: Argument PeakPower must be a real-valued positive scalar! Vectorization is not supported.')
end
if ~isscalar(varianceComplexAWGN) || ~isreal(varianceComplexAWGN) || varianceComplexAWGN < 0
    error('channelCapacity_ComplexAwgnChan_PeakPowerConstr: Argument VarianceOfComplexAWGN must be a real-valued positive scalar! Vectorization is not supported.')
end

SNR_input = peakPower / varianceComplexAWGN;
SNR_input_dB = 10*log10(SNR_input);

% The following load command loads the following variables:
% SNR, SNR_dB, K_evolution, ak_evolution, qk_evolution, R_max, R_LB, R_UD.
% If heavy use of this function is intended, then a different code
% design should be employed, where the file is loaded once globally and
% not once per function call.
load channelCapacity_ComplexAwgnChan_PeakPowerConstr_4to40dB

if SNR_input < 3.011
    rate = computeRate_ComplexAwgnChan_ConcentricCircles(1/SNR_input, 1, 1);
    K = 1;
    ak = sqrt(peakPower);
    qk = 1;
    return
end

if nargout > 1 && SNR_input_dB > 29 && SNR_input_dB <= SNR_dB(end)
    warning('channelCapacity_ComplexAwgnChan_PeakPowerConstr: For SNR > 29 dB, the circle parameter set might be unreliable (input SNR is %.2f dB)', SNR_input_dB)
end

if SNR_input_dB > SNR_dB(end) && SNR_input_dB <= 80
    warning('channelCapacity_ComplexAwgnChan_PeakPowerConstr: The input SNR of %.2f dB is too large for normal computation. Using the tight lower bound R_UD instead!', SNR_input_dB)
    rate = computeRate_ComplexAwgnChan_UniformDistrOverUnitDisk(1/SNR_input);
    rate = max([R_max(end), rate]); % avoid a sudden jump
    K = nan;
    ak = nan(K_evolution(end),1);
    qk = nan(K_evolution(end),1);
    return
end

if SNR_input_dB > 80
    warning('channelCapacity_ComplexAwgnChan_PeakPowerConstr: The input SNR of %.2f dB is way too large for normal computation. Using the tight lower bound log2(1 + SNR/e) instead!', SNR_input_dB)
    rate = log2(1 + SNR_input / exp(1));
    K = nan;
    ak = nan(K_evolution(end),1);
    qk = nan(K_evolution(end),1);
    return
end

[~,idx1] = min(abs(SNR - SNR_input));
SNR1 = SNR(idx1);
K = K_evolution(idx1);

rate1 = R_max(idx1);
ak1 = ak_evolution(1:K,idx1);
qk1 = qk_evolution(1:K,idx1);

if abs(SNR1 - SNR_input) < eps
    % The precomputed data has an extremely close match for the input SNR.
    % So we just return that data without any need for interpolation:
    rate = rate1;
    ak = sqrt(peakPower) * ak1;
    qk = qk1;
    return
end

foundTwoNeighboringDataPoints = true;
naturalIncrement = sign(SNR_input - SNR1);
if K_evolution(idx1 + naturalIncrement) == K
    % That's perfect. Let's go for interpolation!
    idx2 = idx1 + naturalIncrement;
elseif K_evolution(idx1 - naturalIncrement) == K
    % That's okay too. Let's go for extrapolation!
    idx2 = idx1 - naturalIncrement;
else
    warning('channelCapacity_ComplexAwgnChan_PeakPowerConstr: There is only one data point with K = %d circles, at SNR = %.2f dB (input SNR is %.2f dB). Can not interpolate or extrapolate. Returning ak, qk for the single neighboring data point.', K, 10*log10(SNR1), SNR_input_dB)
    idx2 = idx1 + naturalIncrement;
    foundTwoNeighboringDataPoints = false;
    ak = sqrt(peakPower) * ak1;
    qk = qk1;
end

SNR2 = SNR(idx2);
rate2 = R_max(idx2);

slope = (rate2 - rate1) / (SNR2 - SNR1);
rate = slope * (SNR_input - SNR1) + rate1;

if foundTwoNeighboringDataPoints
    ak2 = ak_evolution(1:K,idx2);
    slope = (ak2 - ak1) ./ (SNR2 - SNR1);
    ak = slope * (SNR_input - SNR1) + ak1;
    ak = sqrt(peakPower) * ak;
    
    qk2 = qk_evolution(1:K,idx2);
    slope = (qk2 - qk1) ./ (SNR2 - SNR1);
    qk = slope * (SNR_input - SNR1) + qk1;
    
    % Extrapolation can *sometimes* cause a negative qk(end). Fix that!
    qk(qk <= 0) = eps;
    qk = qk / sum(qk);
end
