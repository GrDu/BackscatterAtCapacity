% Computes the achievable information rates of transmit signaling with
% finite constellations, namely with PSK, QAM, and APSK. It plots these
% rates (and reference rates, e.g., the channel capacity) versus the
% signal-to-noise ratio (SNR).
%
% This generates Fig.6 (with some add-ons) of the paper stated below.
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

SNRTargets_APSK_dB = [9, 15, 21];
M_targets = [16, 64, 256];
SNR_targets_APSK = 10.^(SNRTargets_APSK_dB/10);
styles = {'-.','--',':'};
marks = {'o','d','v'};

SNR_dB_min = -10;
SNR_dB_max = 37;
SNR_dB = SNR_dB_min : .5 : SNR_dB_max;
SNR = 10.^(SNR_dB / 10);

[R_APSK, R_QAM, R_PSK] = deal(nan(length(M_targets),length(SNR)));
figure(60101), clf
for idxM = 1 : length(M_targets)
    M = M_targets(idxM)

    [S_APSK, pmf_APSK] = getConstellationAPSK_OptimizedGivenM(M, SNR_targets_APSK(idxM));
    S_QAM = getConstellationQAM(M_targets(idxM));
    S_PSK = getConstellationPSK(M_targets(idxM));

    H_source_APSK = computeEntropyFromPMF(pmf_APSK)
    
    % Plot the symbol constellations
    subplot(3,3,idxM), plot(real(S_APSK),imag(S_APSK),'k.','LineWidth',2)
    axis equal, grid on, xlim([-1.1,1.1]), ylim([-1.1,1.1])
    subplot(3,3,idxM+3), plot(real(S_QAM),imag(S_QAM),'r.','LineWidth',2)
    axis equal, grid on, xlim([-1.1,1.1]), ylim([-1.1,1.1])
    subplot(3,3,idxM+6), plot(real(S_PSK),imag(S_PSK),'b.','LineWidth',2)
    axis equal, grid on, xlim([-1.1,1.1]), ylim([-1.1,1.1])
    
    for n = 1 : length(SNR)
        R_Max = channelCapacity_ComplexAwgnChan_PeakPowerConstr(1, 1/SNR(n));
        R_MaxPSK = computeRate_ComplexAwgnChan_ConcentricCircles(1/SNR(n), 1, 1);
        R_APSK(idxM,n) = computeRate_ComplexAwgnChan_MassPoints(1/SNR(n), S_APSK, pmf_APSK);
        R_QAM (idxM,n) = computeRate_ComplexAwgnChan_MassPoints(1/SNR(n), S_QAM);
        R_PSK (idxM,n) = computeRate_ComplexAwgnChan_MassPoints(1/SNR(n), S_PSK);
    end
end
[R_Max, R_MaxPSK] = deal(nan(1,length(SNR)));
for n = 1 : length(SNR)
    R_Max(n) = channelCapacity_ComplexAwgnChan_PeakPowerConstr(1, 1/SNR(n));
    R_MaxPSK(n) = computeRate_ComplexAwgnChan_ConcentricCircles(1/SNR(n), 1, 1);
end

% Plot all rates vs SNR in the same plot (it's gonna be messy)
figure(6), clf
plot(SNR_dB,R_Max,'b','LineWidth',1,'DisplayName','channel capacity, general passive load'), hold on
plot(SNR_dB,R_MaxPSK,'k','LineWidth',1,'DisplayName','channel capacity, purely reactive load')
for idxM = 1 : length(M_targets)
    plot(SNR_dB,R_PSK(idxM,:),'k','LineStyle',styles{idxM},'LineWidth',1.5,'DisplayName',sprintf('%d-PSK',M_targets(idxM)))
end
for idxM = 1 : length(M_targets)
    plot(SNR_dB,R_QAM(idxM,:),'r','LineStyle',styles{idxM},'LineWidth',1.5,'DisplayName',sprintf('%d-QAM',M_targets(idxM)))
end
for idxM = 1 : length(M_targets)
    plot(SNR_dB,R_APSK(idxM,:),'Color',[0,.8,0],'LineStyle',styles{idxM},'LineWidth',1.5,'HandleVisibility','off')
    [~,nTarget] = min(abs(SNR - SNR_targets_APSK(idxM)));
    plot(SNR_dB(nTarget),R_APSK(idxM,nTarget),'Marker',marks{idxM},'Color',[0,.8,0],'LineWidth',1.5,'HandleVisibility','off')
    plot(nan,nan,'Color',[0,.8,0],'LineStyle',styles{idxM},'Marker',marks{idxM},'LineWidth',1.5,'DisplayName',...
        sprintf('%d-APSK optimized for %d dB',M_targets(idxM),SNRTargets_APSK_dB(idxM)))
end
xlabel('SNR [dB]'), ylabel('achievable information rate [bpcu]')
legend('show', 'Location', 'NorthWest'), grid on
xlim([SNR_dB_min, SNR_dB_max]), ylim([0 10])

% Plot all rates vs SNR again, but nicely distributed across three subplots
figure(60102), clf
subplot(1,3,1)
semilogy(SNR_dB,R_MaxPSK,'k','LineWidth',1,'DisplayName','channel capacity, purely reactive load'), hold on
for idxM = 1 : length(M_targets)
    semilogy(SNR_dB,R_PSK(idxM,:),'k','LineStyle',styles{idxM},'LineWidth',1.5,'DisplayName',sprintf('%d-PSK',M_targets(idxM)))
end
xlabel('SNR [dB]'), ylabel('achievable information rate [bpcu]')
legend('show', 'Location', 'SouthEast'), grid on
xlim([SNR_dB_min, SNR_dB_max]), ylim([1E-1, 1E1])

subplot(1,3,2)
semilogy(SNR_dB,R_Max,'b','LineWidth',1,'DisplayName','channel capacity, general passive load'), hold on
for idxM = 1 : length(M_targets)
    semilogy(SNR_dB,R_QAM(idxM,:),'r','LineStyle',styles{idxM},'LineWidth',1.5,'DisplayName',sprintf('%d-QAM',M_targets(idxM)))
end
xlabel('SNR [dB]'), ylabel('achievable information rate [bpcu]')
legend('show', 'Location', 'SouthEast'), grid on
xlim([SNR_dB_min, SNR_dB_max]), ylim([1E-1, 1E1])

subplot(1,3,3)
semilogy(SNR_dB,R_Max,'b','LineWidth',1,'DisplayName','channel capacity, general passive load'), hold on
for idxM = 1 : length(M_targets)
    plot(SNR_dB,R_APSK(idxM,:),'Color',[0,.8,0],'LineStyle',styles{idxM},'LineWidth',1.5,'HandleVisibility','off')
    [~,nTarget] = min(abs(SNR - SNR_targets_APSK(idxM)));
    plot(SNR_dB(nTarget),R_APSK(idxM,nTarget),'Marker',marks{idxM},'Color',[0,.8,0],'LineWidth',1.5,'HandleVisibility','off')
    plot(nan,nan,'Color',[0,.8,0],'LineStyle',styles{idxM},'Marker',marks{idxM},'LineWidth',1.5,'DisplayName',...
        sprintf('%d-APSK optimized for %d dB',M_targets(idxM),SNRTargets_APSK_dB(idxM)))
end
xlabel('SNR [dB]'), ylabel('achievable information rate [bpcu]')
legend('show', 'Location', 'SouthEast'), grid on
xlim([SNR_dB_min, SNR_dB_max]), ylim([1E-1, 1E1])
