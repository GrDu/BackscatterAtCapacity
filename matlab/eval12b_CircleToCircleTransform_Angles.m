% An evaluation concerning the Möbius transform that underlies the
% Smith chart, mapping circles to other circles.
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

a = .95 : -.15 : .05; % chosen circle radius (s-domain)
%a = [a, 1E-3];
colors = 'rckygbmrckygbm';

NRes = 360;
theta  = linspace(0, 2*pi, NRes+1); % for circles
theta_ = linspace(0, 2*pi, 17); theta_ = theta_(1 : end-1); % for dots

figure(120001), clf

for n = 1 : length(a)
    Gamma  = a(n) * exp(1i*theta);
    Gamma_ = a(n) * exp(1i*theta_);
    
    zL  = (1 + Gamma ) ./ (1 - Gamma );
    zL_ = (1 + Gamma_) ./ (1 - Gamma_);
    zL (abs(zL ) > 100) = nan;
    zL_(abs(zL_) > 100) = nan;

    subplot(2,2,1)
    plot(real(Gamma ), imag(Gamma ), colors(n), 'LineStyle', '-'), hold on
    plot(real(Gamma_), imag(Gamma_), colors(n), 'Marker', '.', 'Markersize', 8, 'LineStyle', 'none')
    grid on, axis equal
    title('circles in reflection coefficient domain')
    xlabel('Re \Gamma')
    ylabel('Im \Gamma')
    xlim([-1.2 1.2])
    ylim([-1.2 1.2])
 
    subplot(2,2,2)
    plot(real(zL ), imag(zL ), colors(n), 'LineStyle', '-'), hold on
    plot(real(zL_), imag(zL_), colors(n), 'Marker', '.', 'Markersize', 8, 'LineStyle', 'none')
    
    grid on,  axis equal
    title('correspondinng circles in impedance domain')
    xlabel('normalized resistance r')
    ylabel('normalized reactance x')
    c = 3;
    xlim([0 2*c]), ylim([-c c])
    set(gca, 'XTick',  0 : 2*c)
    set(gca, 'YTick', -c : c)
end

for n = 1 : length(a)
    Gamma = a(n) * exp(1i*theta);
    zL  = (1 + Gamma ) ./ (1 - Gamma );
    
    r0 = (1 + a(n)^2) / (1 - a(n)^2); % circle center in zL-domain
    beta = angle(zL - r0); % makes sense!
    beta(beta < 0) = beta(beta < 0) + 2*pi;
    
    subplot(2,2,3)
    plot(theta / pi, beta / pi, colors(n), 'LineStyle', '-'), hold on
    grid on,  axis equal
    title('angle-to-angle map, by simulation')
    xlabel('\theta / \pi')
    ylabel('\beta / \pi')
    %xlim([min(theta), max(theta)] / pi)
    xlim([0 2])
    ylim([0 2])
    
    % % % equation (4.16) in Johnny SA v1 (original): (CUT b/c s-DOMAIN STUFF)
    
    % from Gregor's alternative definitions
    beta = 2 * atan(sin(theta) ./ (cos(theta) - a(n))) - theta;
    while any(beta < 0), beta(beta < 0) = beta(beta < 0) + 2*pi; end
    
    subplot(2,2,4)
    plot(theta / pi, beta / pi, colors(n), 'LineStyle', '-'), hold on
    grid on,  axis equal
    title('angle-to-angle map, by derived formula')
    xlabel('\theta / \pi')
    ylabel('\beta / \pi')
    %xlim([min(theta), max(theta)] / pi)
    xlim([0 2])
    ylim([0 2])
end

% attempt to compute PDF of beta
figure(120002), clf

for n = 1 : length(a)
    beta = 2 * atan(sin(theta) ./ (cos(theta) - a(n))) - theta;
    while any(beta < 0), beta(beta < 0) = beta(beta < 0) + 2*pi; end
    
    f_theta = (2*pi*(1 - a(n)^2))^-1 * (a(n)^2 - 2*a(n)*cos(theta) + 1);
    
    d_phi_z = beta(2:end) - beta(1:end-1); % "infinitesimal" step size
    approx_sum = sum(d_phi_z .* f_theta(2:end)) % should sum to 1
    
    plot(beta / pi, f_theta * pi, colors(n), 'LineStyle', '-'), hold on
    grid on
    xlabel('\beta / \pi')
    ylabel('PDF \cdot \pi')
    xlim([0 2])
    ylim([0 5])
end

