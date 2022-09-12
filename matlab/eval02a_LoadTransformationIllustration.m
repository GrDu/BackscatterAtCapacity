% Creates a colorful illustration (very similar to the Smith chart) of the
% mapping between the reflection coefficient Gamma and the
% normalized impedance z of a passive one-port network.
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

colors1 = [[.6 .6  1]; ...
           [.2 .2  1]; ...
           [ 0  0 .6]; ...
           [ 0  0  0]; ...
           [.6  0  0]; ...
           [ 1 .2 .2]; ...
           [ 1 .6 .6]];
colors2 = [[ 0  0  0]; ...
           [ 0 .7  0]; ...
           [ 0 .9  0]; ...
           [.5  1 .5]; ...
           [.8  1 .8]];
      
markers = 'o^*.';
markSize = 8;
lw = 1.5;

rVal =  0 : 4; Nr = length(rVal);
xVal = -3 : 3; Nx = length(xVal);

[r,x] = meshgrid(rVal, xVal);
z = r + 1i * x;

Gamma = (z - 1) ./ (z + 1);

figure(20001), clf

subplot(1,2,1)

for idx = 1 : length(rVal)
    rad = 1 ./ (1 + rVal(idx));
    phi = linspace(-pi,pi,61);
    plot(1-rad + rad*cos(phi), rad*sin(phi), 'Color', colors2(idx,:), 'LineWidth', lw, 'HandleVisibility', 'off'), hold on
end
for idx = 1 : length(xVal)
    r_ = [0 : .2 : 6, 10, 20, 50, 100, 1E12];
    z_ = r_ + 1i * xVal(idx);
    Gamma_ = (z_ - 1) ./ (z_ + 1);
    plot(real(Gamma_), imag(Gamma_), '--', 'Color', colors1(idx,:), 'LineWidth', lw, 'HandleVisibility', 'off')
    
    for add = 0 : 3
        Gamma_ = (add + 1i * xVal(idx) - 1) ./ (add + 1i * xVal(idx) + 1);
        plot(real(Gamma_), imag(Gamma_), 'Marker', markers(add+1), 'Color', colors1(idx,:), 'MarkerSize', markSize, 'LineWidth', lw, 'HandleVisibility', 'off')
    end
end
axis equal
set(gca, 'Box', 'off')
xlim([-1 1])
ylim([-1 1])
xlabel('Re \Gamma')
ylabel('Im \Gamma')

subplot(1,2,2)

for idx = 1 : length(rVal)
    plot(rVal(idx)*[1 1], [-13 13], 'Color', colors2(idx,:), 'LineWidth', lw, 'HandleVisibility', 'off'), hold on
end
for idx = 1 : length(xVal)
    plot([min(rVal) 13], xVal(idx)*[1 1], '--',    'Color', colors1(idx,:), 'MarkerSize', markSize, 'LineWidth', lw, 'HandleVisibility', 'off')
    for add = 0 : 3
        plot(rVal(add+1), xVal(idx), 'Marker', markers(add+1), 'Color', colors1(idx,:), 'MarkerSize', markSize, 'LineWidth', lw, 'HandleVisibility', 'off')
    end
end
axis equal
set(gca, 'Box', 'off')
xlim([-.5 max(rVal)+.5])
ylim(3.5*[-1 1])
xlabel('Re(z)')
ylabel('Im(z)')
