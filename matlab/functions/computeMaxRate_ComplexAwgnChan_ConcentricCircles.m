function [R_max, ak, qk] = computeMaxRate_ComplexAwgnChan_ConcentricCircles(VarianceComplexNoise, ak_init, qk_init)
%computeMaxRate_ComplexAwgnChan_ConcentricCircles: Maximizes and returns
% the achievable information rate over a complex-valued AWGN channel whereby
% whereby the transmit distribution is supported on a set of concentric
% circles. Herein the circle radii ak and probabilities qk are optimized.
% The number of circles is given and fixed.
%
%   Usage: [R_max, ak, qk] = ...
%   computeMaxRate_ComplexAwgnChan_ConcentricCircles(VarianceComplexNoise, ak_init, qk_init)
%
%   Inputs:
%   VarianceComplexNoise is the variance of the complex-valued AWGN.
%   ak_init is a vector with the initial circle radii.
%   qk_init is a vector with the initial circle probabilities.
%
%   Outputs:
%   R_max is the resulting information rate in bit per channel use (bpcu).
%       This is the channel capacity if the function was
%       initialized appropriately and with the correct number of circles K.
%   ak_init is a vector with the optimized circle radii.
%   qk_init is a vector with the optimized circle probabilities.

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

    if ~issorted(ak_init, 'descend')
        error('maximizeAchievableRateFromCircles: Elements of vector ak_init (the initial circle radii) must be in ascending order!')
    end
    NCircles = length(ak_init);
    
    if NCircles == 1
        ak = ak_init(1); 
        qk = 1;
        R_max = computeRate_ComplexAwgnChan_ConcentricCircles(VarianceComplexNoise, ak, qk);
        return
    end
    
    % stack the inits in one vector for technical reasons
    x_init = [ak_init(2 : end), qk_init];
    
    % objective to minimize (maximize the rate by minimizing the negative rate)
    negativeRate = @(x) ...
       -computeRate_ComplexAwgnChan_ConcentricCircles(VarianceComplexNoise, [ak_init(1), x(1 : NCircles-1)], x(NCircles : end));
    
    % Now use the Matlab function minimization (constrained), a.k.a. fmincon.
    % For details type:
    % doc fmincon
    % help fmincon
    options = optimoptions('fmincon','Algorithm','interior-point',...
        'Display', 'none', ... % 'final-detailed'
        'OptimalityTolerance', 100*eps);
        %'StepTolerance', 100*eps, ...
        %'MaxIterations', 15E3, ...
        %'MaxFunctionEvaluations', 50E3);
    
    % prepare matrix with inequality constraints (which preserve the radius order)
    NDim = 2*NCircles - 1; % dimensionality of the optimization problem
    numericalMargin = 1 - 1/(5*NCircles); % e.g., value .9 means that a circle must have radius <= .9 * radius_next
    A = diag(ones(1,NCircles-2), 1) - numericalMargin * eye(NCircles-1);
    A = [A(1 : end-1,:), zeros(NCircles-2,NCircles)]; % discard last row
    
    %testCostFunctionValue = negativeRate(x_init) % sanity check debug output
    
    [x_opt, negativeRate_opt, exitflag, output] = fmincon(negativeRate, ...
        x_init,                                     ... % initial value for the iterative function minimization
        A, zeros(NCircles-2,1),                     ... % linear inequality constraints: preserve the radius order
        [zeros(1,NCircles-1), ones(1,NCircles)], 1, ... % linear equality constraint: sum of PMF must be 1
        zeros(1,NDim),                              ... % lower bounds: 0 <= radius, 0 <= probability
        [numericalMargin*ak_init(1)*ones(1,NCircles-1), ones(1,NCircles)], ... % upper bounds: radius <= 1, probability <= 1
        [],                                         ... % don't use any nonlinear constraints
        options);
    
    R_max = -negativeRate_opt;
    ak = [ak_init(1), x_opt(1 : NCircles-1)];
    qk = x_opt(NCircles : end);
    
    % % diagnostic output
    % output, exitflag, C, radii_opt, pmf_opt
end
