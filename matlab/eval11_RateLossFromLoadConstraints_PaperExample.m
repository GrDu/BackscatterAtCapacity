% Computes the high-SNR rate loss of technical constraints on a backscatter
% tag load and generates a useful illustration.
%
% This generates Fig.11 of the paper stated below.
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

%% Configure the experiment

Delta = 0.25;
Q = 10;

%% Helper quantities (see paper)

x_lb = -Q*Delta/(1-Delta);
x_ub =  Q*Delta/(1+Delta);
rVal = [0, logspace(-3, 5, 41)]; % load resistance values: a wide range that
                                 % is representative of the positive real numbers

%% Prepare the figure

myFig = figure(11), clf
axis equal, hold on
xlabel('Re \Gamma'), ylabel('Im \Gamma')
xlim([-1 1]), ylim([-1 1])

%% Investigate the part above the load reactance upper constraint

xVal = sort([linspace(x_ub, 8*x_ub, 61), logspace(log10(x_ub)+1, log10(x_ub)+9, 33)]);
z1 = 1i * xVal;
z2 = rVal + 1i * max(xVal);
z3 = max(rVal) + 1i * xVal(end : -1 : 1);
z4 = rVal(end : -1 : 1) + 1i * min(xVal);
z = [z1, z2(2:end), z3(2:end), z4(2:end)];
Gamma = (z-1) ./ (z+1);

%% Plot the result

patch(real(Gamma), imag(Gamma),'r')
plot(real(Gamma), imag(Gamma),'k-','LineWidth',2)

%% Investigate the part below the load reactance lower constraint

xVal = sort([linspace(8*x_lb,x_lb,61), -logspace(log10(abs(x_lb))+1, log10(abs(x_lb))+9,33)]);
z1 = 1i * xVal;
z2 = rVal + 1i * max(xVal);
z3 = max(rVal) + 1i * xVal(end : -1 : 1);
z4 = rVal(end : -1 : 1) + 1i * min(xVal);
z = [z1, z2(2:end), z3(2:end), z4(2:end)];
Gamma = (z-1) ./ (z+1);

%% Plot the result

patch(real(Gamma), imag(Gamma),'r')
plot(real(Gamma), imag(Gamma),'k-','LineWidth',2)

%% Plot the bounding circle 

phi = linspace(-pi,pi,121);
plot(cos(phi), sin(phi), 'k-', 'LineWidth', 2)

%% Compute the relative constraint-restricted area and the resulting high-SNR rate loss

countAll = 0;
countOut_ub = 0;
countOut_lb = 0;
for GammaRe = -1 : .001 : 1
    for GammaIm = -1 : .001 : 1
        Gamma = GammaRe + 1i * GammaIm;
        if abs(Gamma) > 1
            continue;
        end
        
        countAll = countAll + 1;
        z = (1 + Gamma) / (1 - Gamma);
        
        if imag(z) > x_ub
            countOut_ub = countOut_ub + 1;
        elseif imag(z) < x_lb
            countOut_lb = countOut_lb + 1;
        end
    end
end

loss_relative_area = (countOut_ub + countOut_lb) / countAll
loss_bpcu_highSNR = log2(1-loss_relative_area)
