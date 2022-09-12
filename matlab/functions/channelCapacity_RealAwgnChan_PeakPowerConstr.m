function [rate, MOnesided, smOnesided, pmOnesided] = ...
   channelCapacity_RealAwgnChan_PeakPowerConstr(peakPower, varianceRealAWGN)
%channelCapacity_RealAwgnChan_PeakPowerConstr: Returns the channel capacity
% over a real-valued AWGN channel with a constraint on the signal peak-power,
% but no constraint on the average transmit power. It can also return the
% parameters of the (finite number of) mass points that describe the
% capacity-achieving transmit distribution.
%
%   Usage: Obtain just the rate with
%   rate = channelCapacity_RealAwgnChan_PeakPowerConstr(peakPower, varianceRealAWGN)
%   or the rate and mass-point parameters with
%   [rate, MOnesided, smOnesided, pmOnesided] = ...
%      channelCapacity_RealAwgnChan_PeakPowerConstr(peakPower, varianceRealAWGN)
%   
%   Inputs:
%   PeakPower is the maximum peak power in normal scale (not in dB).
%   varianceRealAWGN is the variance of the real-valued additive white 
%       Gaussian noise (AWGN).
%
%   Outputs:
%   rate is the channel capacity in bit per channel use (bpcu), cf. the
%       paper Sec. III-E (Resistive Load Modulation).
%   MOnesided is the number of mass points with location >= 0. These 
%       describe the capacity-achieving transmit distribution.
%   smOnesided is a length-MOnesided vector holding the positive mass point
%       locations. Their negatives are also mass points.
%   pmOnesided is a length-MOnesided vector holding the mass point
%       probabilities. These have to be split in half externally, between
%       positive and negative mass point location.
%
%   This function relies on loading and interpolating precomputed data.
%   Computing everything from scratch would be quite an ordeal because it
%   requires a fine-grained iteration from low SNR to the specified SNR.
%   This process is implemented in the script
%   eval04c_RealAwgnChan_OptMassPointParamsVsSNR_Compute.m
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
      error('channelCapacity_RealAwgnChan_PeakPowerConstr: Argument PeakPower must be a real-valued positive scalar! Vectorization is not supported.')
   end
   if ~isscalar(varianceRealAWGN) || ~isreal(varianceRealAWGN) || varianceRealAWGN < 0
      error('channelCapacity_RealAwgnChan_PeakPowerConstr: Argument VarianceOfComplexAWGN must be a real-valued positive scalar! Vectorization is not supported.')
   end
   
   SNR_input = peakPower / varianceRealAWGN;
   SNR_input_dB = 10*log10(SNR_input);

   % The following load command loads the following variables:
   % SNR_real, SNR_real_dB, MOnesided_evolution, smOnesided_evolution, 
   % pmOnesided_evolution, varianceRealNoise, R_max_real, R_LB_real, R_UD_real.
   % If heavy use of this function is intended, then a different code
   % design should be employed, where the file is loaded once globally and
   % not once per function call.
   load channelCapacity_RealAwgnChan_PeakPowerConstr_0to40dB
   
   if SNR_input_dB < SNR_real_dB(1)
      rate = computeRate_RealAwgnChan_MassPoints(1/SNR_input, [-1 1], [.5 .5]);
      MOnesided = 1;
      smOnesided = sqrt(peakPower);
      pmOnesided = 1;
      return
   end

   if SNR_input_dB > SNR_real_dB(end) && SNR_input_dB <= 80
      warning('channelCapacity_RealAwgnChan_PeakPowerConstr: The input SNR of %.2f dB is too large for normal computation. Using the tight lower bound R_UD instead!', SNR_input_dB)
      rate = computeRate_RealAwgnChan_UniformDistrOverUnitInterval(1/SNR_input);
      rate = max([R_max_real(end), rate]); % avoid a sudden jump
      MOnesided = nan;
      smOnesided = nan(MOnesided_evolution(end),1);
      pmOnesided = nan(MOnesided_evolution(end),1);
      return
   end
   
   if SNR_input_dB > 80
      warning('channelCapacity_RealAwgnChan_PeakPowerConstr: The input SNR of %.2f dB is way too large for normal computation. Using the tight lower bound log2(1 + SNR/e) instead!', SNR_input_dB)
      rate = .5 * log2(1 + 4*SNR_input/(pi*exp(1)));
      MOnesided = nan;
      smOnesided = nan(MOnesided_evolution(end),1);
      pmOnesided = nan(MOnesided_evolution(end),1);
      return
   end
   
   [~,idx1] = min(abs(SNR_real - SNR_input));
   SNR1 = SNR_real(idx1);
   MOnesided = MOnesided_evolution(idx1);
   
   rate1 = R_max_real(idx1);
   smOnesided1 = smOnesided_evolution(1:MOnesided,idx1);
   pmOnesided1 = pmOnesided_evolution(1:MOnesided,idx1);
   
   if abs(SNR1 - SNR_input) < eps
      % The precomputed data has an extremely close match for the input SNR.
      % So we just return that data without any need for interpolation:
      rate = rate1;
      smOnesided = sqrt(peakPower) * smOnesided1;
      pmOnesided = pmOnesided1;
      return
   end
   
   foundTwoNeighboringDataPoints = true;
   naturalIncrement = sign(SNR_input - SNR1);
   if MOnesided_evolution(idx1 + naturalIncrement) == MOnesided
      % That's perfect. Let's go for interpolation!
      idx2 = idx1 + naturalIncrement;
   elseif MOnesided_evolution(idx1 - naturalIncrement) == MOnesided
      % That's okay too. Let's go for extrapolation!
      idx2 = idx1 - naturalIncrement;
   else
      warning('channelCapacity_RealAwgnChan_PeakPowerConstr: There is only one data point with MOnesided = %d mass points, at SNR = %.2f dB (input SNR is %.2f dB). Can not interpolate or extrapolate. Returning sm, pm for the single neighboring data point.', MOnesided, 10*log10(SNR1), SNR_input_dB)
      idx2 = idx1 + naturalIncrement;
      foundTwoNeighboringDataPoints = false;
      smOnesided = sqrt(peakPower) * smOnesided1;
      pmOnesided = pmOnesided1;
   end
   
   SNR2 = SNR_real(idx2);
   rate2 = R_max_real(idx2);
   
   slope = (rate2 - rate1) / (SNR2 - SNR1);
   rate = slope * (SNR_input - SNR1) + rate1;
   
   if foundTwoNeighboringDataPoints
      ak2 = smOnesided_evolution(1:MOnesided,idx2);
      slope = (ak2 - smOnesided1) ./ (SNR2 - SNR1);
      smOnesided = slope * (SNR_input - SNR1) + smOnesided1;
      smOnesided = sqrt(peakPower) * smOnesided;
   
      qk2 = pmOnesided_evolution(1:MOnesided,idx2);
      slope = (qk2 - pmOnesided1) ./ (SNR2 - SNR1);
      pmOnesided = slope * (SNR_input - SNR1) + pmOnesided1;
   end
end
