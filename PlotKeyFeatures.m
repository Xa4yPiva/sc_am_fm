function [] = PlotKeyFeatures(kf, ci, snr, modNames)

% figure('NumberTitle', 'off', 'Name', '4 params');
% set(gca, 'DefaultAxesFontSize', 18); set(gcf, 'color', 'w');
% subplot(2,2,1);
% title('\gamma_{max}');
% xlabel('SNR, dB');
% for i = 1 : size(kf, 1)
%     hold on; plot(snr, [kf(i,:).gammaMax], 'linewidth', 2);
% end
% grid on;
% legend(modNames, 'location', 'northwest'); legend('show');
% 
% subplot(2,2,2);
% title('\sigma_{ap}');
% xlabel('SNR, dB');
% for i = 1 : size(kf, 1)
%     hold on; plot(snr, [kf(i,:).sigmaAP], 'linewidth', 2);
% end
% grid on;
% legend(modNames, 'location', 'southwest'); legend('show');
% 
% subplot(2,2,3);
% title('\sigma_{dp}');
% xlabel('SNR, dB');
% for i = 1 : size(kf, 1)
%     hold on; plot(snr, [kf(i,:).sigmaDP], 'linewidth', 2);
% end
% grid on;
% legend(modNames, 'location', 'southwest'); legend('show');
% 
% subplot(2,2,4);
% title('P');
% xlabel('SNR, dB');
% for i = 1 : size(kf, 1)
%     hold on; plot(snr, [kf(i,:).P], 'linewidth', 2);
% end
% grid on;
% legend(modNames, 'location', 'northwest'); legend('show');

figure('NumberTitle', 'off', 'Name', 'GammaMax');
set(gca, 'DefaultAxesFontSize', 18); set(gcf, 'color', 'w');
title('\gamma_{max}');
xlabel('SNR, dB');
for i = 1 : size(kf, 1)
    hold on; 
    err = [ci(i,:).gammaMax];
    gammaMax = [kf(i,:).gammaMax];
    errorbar(snr, gammaMax, err(1,:), err(2,:), 'linewidth', 2);
end
grid on;
legend(modNames, 'location', 'northwest'); legend('show');

% figure('NumberTitle', 'off', 'Name', 'SigmaAP');
% set(gca, 'DefaultAxesFontSize', 18); set(gcf, 'color', 'w');
% title('\sigma_{ap}');
% xlabel('SNR, dB');
% for i = 1 : size(kf, 1)
%     hold on; 
%     plot(snr, [kf(i,:).sigmaAP], 'linewidth', 2);
% end
% grid on;
% legend(modNames, 'location', 'southwest'); legend('show');
% 
% figure('NumberTitle', 'off', 'Name', 'SigmaAF');
% set(gca, 'DefaultAxesFontSize', 18); set(gcf, 'color', 'w');
% title('\sigma_{dp}');
% xlabel('SNR, dB');
% for i = 1 : size(kf, 1)
%     hold on; 
%     plot(snr, [kf(i,:).sigmaDP], 'linewidth', 2);
% end
% grid on;
% legend(modNames, 'location', 'southwest'); legend('show');
% 
% figure('NumberTitle', 'off', 'Name', 'P');
% set(gca, 'DefaultAxesFontSize', 18); set(gcf, 'color', 'w');
% title('P');
% xlabel('SNR, dB');
% for i = 1 : size(kf, 1)
%     hold on; 
%     plot(snr, [kf(i,:).P], 'linewidth', 2);
% end
% grid on;
% legend(modNames, 'location', 'southwest'); legend('show');
% 
% figure('NumberTitle', 'off', 'Name', 'SigmaAF');
% set(gca, 'DefaultAxesFontSize', 18); set(gcf, 'color', 'w');
% title('\sigma_{af}');
% xlabel('SNR, dB');
% for i = 1 : size(kf, 1)
%     hold on; 
%     plot(snr, [kf(i,:).sigmaAF], 'linewidth', 2);
% end
% grid on;
% legend(modNames, 'location', 'southwest'); legend('show');
% 
% figure('NumberTitle', 'off', 'Name', 'SigmaDF');
% set(gca, 'DefaultAxesFontSize', 18); set(gcf, 'color', 'w');
% title('\sigma_{df}');
% xlabel('SNR, dB');
% for i = 1 : size(kf, 1)
%     hold on; 
%     plot(snr, [kf(i,:).sigmaDF], 'linewidth', 2);
% end
% grid on;
% legend(modNames, 'location', 'southwest'); legend('show');

end

