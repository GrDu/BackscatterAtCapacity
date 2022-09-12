function [R_max_real, smOnesided, pmOnesided] = computeMaxRate_RealAwgnChan_MassPoints(varianceRealNoise, smOnesided_init, pmOnesided_init)
%computeMaxRate_RealAwgnChan_MassPoints: Optimizes the parameters of
% the "positive-side" mass points that describe a symmetric discrete
% transmit distribution, to maximize the information rate over
% a real-valued peak-power-limited AWGN channel .
%
%   Usage: [R_max_real, smOnesided, pmOnesided] = ...
%   computeMaxRate_RealAwgnChan_MassPoints(varianceRealNoise, smOnesided_init, pmOnesided_init)
%
%   Inputs:
%   varianceRealNoise is the AWGN variance.
%   smOnesided_init is a vector with the initial locations
%       of those mass points with positive location.
%   pmOnesided_init is a vector with the initial probabilities
%       of those mass points with positive location.
%   Outputs:
%   R_max_real is the resulting information rate in bit per channel use (bpcu).
%   smOnesided is a vector with the optimized locations of those
%       mass points with positive location.
%   pmOnesided is a vector with the optimized probabilities of those
%       mass points with positive location.
%
%   All entries of smOnesided_init, smOnesided must/will be between 0 and 1.
%   In other words, the support of the transmit distribution is assumed to
%   be a subset of the interval [-1,1].
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

    if ~issorted(smOnesided_init, 'descend')
        error('computeMaxRate_RealAwgnChan_MassPoints: Elements of vector radii_init must be in ascending order!')
    end
    MOnesided = length(smOnesided_init);
    
    if MOnesided == 1
        smOnesided = smOnesided_init(1);
        pmOnesided = 1;
        R_max_real = computeRateGivenOnesidedData(varianceRealNoise, smOnesided, pmOnesided);
        return
    end
    
    % stack the inits in one vector for technical reasons
    x_init = [smOnesided_init(2 : end), pmOnesided_init];
    
    % objective to minimize (maximize the rate by minimizing the negative rate)
    negativeRate = @(x) -computeRateGivenOnesidedData(varianceRealNoise, [smOnesided_init(1), x(1 : MOnesided-1)], x(MOnesided : end));
    
    % Now use the Matlab function minimization (constrained), a.k.a. fmincon.
    % For details type:
    % doc fmincon
    % help fmincon
    options = optimoptions('fmincon','Algorithm','interior-point',...
        'Display', 'none'); ... % 'final-detailed'
        %'OptimalityTolerance', 100*eps);
        %'StepTolerance', 100*eps, ...
        %'MaxIterations', 15E3, ...
        %'MaxFunctionEvaluations', 50E3);
    
    % prepare inequality constraints matrix (preserve radius order)
    NDim = 2*MOnesided - 1; % dimensionality of optimization problem
    numericalMargin = 1 - 1/(5*MOnesided); % e.g., value .9 means that a positive point must have sm <= .9 * sm_next
    A = diag(ones(1,MOnesided-2), 1) - numericalMargin * eye(MOnesided-1);
    A = [A(1 : end-1,:), zeros(MOnesided-2,MOnesided)]; % discard last row
    
    [x_opt, negativeRate_opt, exitflag, output] = fmincon(negativeRate, ...
        x_init,                                     ... % initial value for the iterative function minimization
        A, zeros(MOnesided-2,1),                     ... % linear inequality constraints: preserve smOnesided order
        [zeros(1,MOnesided-1), ones(1,MOnesided)], 1, ... % linear equality constraint: sum of PMF must be 1
        zeros(1,NDim),                              ... % lower bounds: 0 <= smOnesided, 0 <= pm
        [numericalMargin*smOnesided_init(1)*ones(1,MOnesided-1), ones(1,MOnesided)], ... % upper bounds: smOnesided <= 1, pm <= 1
        [],                                         ... % don't use any nonlinear constraints
        options);
    
    R_max_real = -negativeRate_opt;
    smOnesided = [smOnesided_init(1), x_opt(1 : MOnesided-1)];
    pmOnesided = x_opt(MOnesided : end);
    
    % % diagnostic output
    % output, exitflag, C, smOnesided, pmOnesided
end

function rate = computeRateGivenOnesidedData(VarOfRealNoise, smOnesided, pmOnesided)
   if any(smOnesided < 0)
      error('computeRateGivenOnesidedData: expects only positive values input vector smOnesided (mirroring is done internally)')
   end
   
   % Do the mirroring:
   sm = [-smOnesided(end:-1:1), smOnesided];
   pm = .5 * [pmOnesided(end:-1:1), pmOnesided];
   
   % Compute achievable rate of constellation:
   rate = computeRate_RealAwgnChan_MassPoints(VarOfRealNoise, sm, pm);
end
