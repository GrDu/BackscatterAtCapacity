% Generates the intensity plots of the probability density functions
% (in terms of real and imaginary parts) of the complex-valued normalized
% impedance of a backscatter tag load. This distribution guarantees that
% the "transmit signal", i.e. the load reflection coefficient, has maximum
% entropy. This yields near-capacity rates at high SNR.
%
% This generates Fig.13 of the paper stated below.
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

NrVal = 402; % Avoid to hit z=1 precisely to avoid pdf=NaN there
rVal = linspace(-1, 3, NrVal);
NxVal = 401;
xVal = linspace(-2, 2, NxVal);

[r,x] = meshgrid(rVal,xVal);
z = r + 1i * x;
Gamma = (z - 1) ./ (z + 1);
a_  = abs(Gamma);
th_ = angle(Gamma);

pdf = nan(size(z));
for nrVal = 1 : NrVal
    for nxVal = 1 : NxVal
        if rVal(nrVal) < 0
            pdf(nxVal,nrVal) = 0;
            continue
        end
        
        a  = a_(nxVal,nrVal);
        th = th_(nxVal,nrVal);
        
        N = 1 + a^2 - 2 * a * cos(th);
        dNda  = 2 * (a - cos(th));
        dNdth = 2 * a * sin(th);
        M_tmp = 2 * N * [-a, 0; sin(th), a * cos(th)];
        v_tmp = [1-a^2;  2*a*sin(th)];
        J = 1/N^2 * (M_tmp - v_tmp * [dNda, dNdth]);
        
        pdf(nxVal,nrVal) = a/pi / abs(det(J));
    end
end

figure(13), clf
imagesc(rVal,xVal,pdf)
%contour(rVal,xVal,pdf, .1 : .3 : 1)
colormap(flipud(gray))
hold on
contour(rVal,xVal,pdf)
hold off
grid on
axis equal
xlim([min(rVal), max(rVal)]), xlabel('normalized load resistance r'), xticks(min(rVal) : max(rVal))
ylim([min(xVal), max(xVal)]), ylabel('normalized load reactance x'),  yticks(min(xVal) : max(xVal))
h = colorbar;
set(h, 'Ticks', [0 max(pdf(:))])
set(h, 'TickLabels', {'0','4/\pi'})
h = ylabel(h, 'probability density');

% [Sanity Check] How large is the double integral of the PDF over the chosen
% canvas? This shows the heavy tails. For a large canvas, the result should approach 1. 
dr = (max(rVal) - min(rVal)) / (length(rVal)-1);
dx = (max(xVal) - min(xVal)) / (length(xVal)-1);
probabilityMassInCanvas = sum(pdf(:)) * dr * dx

% % plot the PDF along the real axis
% figure(130001), clf
% plot(rVal,pdf((NxVal+1)/2,:))
