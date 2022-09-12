% Demonstrates a simple known construction of a
% amplitude-and-phase-shift keying (APSK) constellation from a
% quadrature-amplitude-modulation (QAM) constellation.
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

M = 64; % choose from 16, 32, 64, 128, 256

% Standard QAM
S_QAM = getConstellationQAM(M);

% According to paper "QAM to Circular Isomorphic Constellations"
% by Kayhan, IEEE 2016
S_ToCirc = nan(M,1);
for m = 1 : M
    s = S_QAM(m);
    S_ToCirc(m) = max(abs([real(s),imag(s)])) * sqrt(2) * s / abs(s);
end

% APSK acoording to DVB-S2X standard,
% https://ch.mathworks.com/help/comm/ref/dvbsapskmod.html
S_DVBS2X = dvbsapskmod(0 : M-1, M, 's2x');

figure(70201), clf
phi = linspace(-pi,pi,181);

subplot(1,3,1)
plot(cos(phi),sin(phi),'k-'), hold on
plot(real(S_QAM),imag(S_QAM),'rx','LineWidth',2,'MarkerSize',8)
title('QAM')
axis equal, grid on
xlim([-1.2,1.2])
ylim([-1.2,1.2])
subplot(1,3,2)
plot(cos(phi),sin(phi),'k-'), hold on
plot(real(S_ToCirc),imag(S_ToCirc),'bx','LineWidth',2,'MarkerSize',8)
title('APSK from QAM (Circular Isomorphic Mapping)')
axis equal, grid on
xlim([-1.2,1.2])
ylim([-1.2,1.2])
subplot(1,3,3)
plot(cos(phi),sin(phi),'k-'), hold on
plot(real(S_DVBS2X),imag(S_DVBS2X),'gx','LineWidth',2,'MarkerSize',8)
title('APSK from DVB-S2X Standard')
axis equal, grid on
xlim([-1.2,1.2])
ylim([-1.2,1.2])
