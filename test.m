addpath(genpath('../sc_common'));
addpath(genpath('../matlab_utils'));
clear;
close all;

%% Signal
fs0 = 12e3;
T = 1;
fLowHz = 300;
fHighHz = 3600;
bandHz = fHighHz - fLowHz;
data = RandomBandLimitedSignal(fs0, T, 20, fLowHz, fHighHz, 4000, 60, 1, 60, 'uniform');
fs = 200e3;
factor = fs / fs0;
[p, q] = rat(factor);
x = resample(data, p, q);
lenSignal = length(x);
N = 2 ^ nextpow2(lenSignal);
freqs = (-N/2 : (N-1)/2) * fs / N;
fc = 60e3;
offset = exp(1i * 2*pi*fc * (0:lenSignal-1)/fs);

%% Modulation
mAM1 = 0.3;
xAM1 = ammod(mAM1 * x, fc, fs, 0, 0.5);
mAM2 = 0.5;
xAM2 = ammod(mAM2 * x, fc, fs, 0, 0.5);
mAM3 = 0.7;
xAM3 = ammod(mAM3 * x, fc, fs, 0, 0.5);
specAM1 = fftshift(abs(fft(xAM1, N))) / (N/2);
xDSB = ammod(x, fc, fs, 0, 0);
specDSB = fftshift(abs(fft(xDSB, N))) / (N/2);
xLSB = ssbmod(x, fc, fs, 0);
specLSB = fftshift(abs(fft(xLSB, N))) / (N/2);
xUSB = ssbmod(x, fc, fs, 0, 'upper');
specUSB = fftshift(abs(fft(xUSB, N))) / (N/2);
fDev1 = 5e3;
xFM1 = fmmod(x, fc, fs, fDev1);
% specFM = fftshift(abs(fft(xFM, N))) / (N/2);
fDev2 = 10e3;
xFM2 = fmmod(x, fc, fs, fDev2);
fDev3 = 25e3;
xFM3 = fmmod(x, fc, fs, fDev3);

xNoise = zeros(1, lenSignal);

% plot(xAM1(1:1024));
% return;
% figure(1);
% subplot(5,1,1); plot(freqs, specAM); grid on;
% subplot(5,1,2); plot(freqs, specDSB); grid on;
% subplot(5,1,3); plot(freqs, specLSB); grid on;
% subplot(5,1,4); plot(freqs, specUSB); grid on;
% subplot(5,1,5); plot(freqs, specFM); grid on;

signals = [xAM1; xDSB; xLSB; xUSB; xFM3;];
% signals = [xAM1; xAM2; xAM3; xDSB; xLSB; xUSB; xFM1; xFM2; xFM3;];
sigsNum = size(signals, 1);
modNames = ["AM-0.3", "DSB", "LSB", "USB", "FM - 25kHz", "Noise"];
% modNames = ["AM - 0.3", "AM - 0.5", "AM - 0.7", ...
%     "DSB", "LSB", "USB", ...
%     strcat("FM - ", num2str(fDev1/1e3), "kHz"), ...
%     strcat("FM - ", num2str(fDev2/1e3), "kHz"), ...
%     strcat("FM - ", num2str(fDev3/1e3), "kHz"), "Noise"];
envelopes = zeros(size(signals));
offsets = [-fc, -fc, -fc + bandHz/2 + fLowHz, -fc - bandHz/2 - fLowHz, -fc, 0];
% offsets = [-fc * ones(1, sigsNum-1), 0];
% offsets = -fc * ones(1, sigsNum);
expOff = exp(1i * 2*pi*offsets' .* (0:lenSignal-1)/fs);
for i = 1 : sigsNum
    envelopes(i, :) = hilbert(signals(i, :)) .* expOff(i, :);
end

fs2 = 80e3;
factor = fs2 / fs;
[p, q] = rat(factor);
envelopes = (resample(envelopes', p, q))';

N = 2 ^ nextpow2(size(envelopes, 2));
freqs = (-N/2 : (N-1)/2) * fs2 / N;
specEnv = fftshift(abs(fft(envelopes, N, 2)) / (N/2), 2);

% figure(2);
% subplot(5,1,1); plot(freqs, mag2db((specEnv(1,:)))); grid on;
% subplot(5,1,2); plot(freqs, mag2db((specEnv(2,:)))); grid on;
% subplot(5,1,3); plot(freqs, mag2db((specEnv(3,:)))); grid on;
% subplot(5,1,4); plot(freqs, mag2db((specEnv(4,:)))); grid on;
% subplot(5,1,5); plot(freqs, mag2db((specEnv(5,:)))); grid on;

%% KeyFeatures
thresholds.a = 1;
lenFrame = 4096;
% kf = KeyFeatures(awgn(envelopes(4, 328 : 328 + lenFrame), 10, 'measured'), thresholds.a)
% return

lenEnv = size(envelopes, 2);
threshA = thresholds.a;
snr = -15 : 1 : 15;
expNum = 100;
iteration = 0;
beta = 0.9;
% h = waitbar(0, 'Computing...');
% cyclesNum = length(snr) * sigsNum * expNum;
disp('Computing Key Features vs SNR ...');
tic
for j = 1 : length(snr)
    snrCurr = snr(j);
    for i = 1 : sigsNum
        env = envelopes(i,:);
%         lenFrame = 2 .^ (7 : 15);
%         for l = 1 : length(lenFrame)
%             for k = 1 : expNum
%                 kf_d(k, l) = KeyFeatures(awgn(env(1:lenFrame(l)), snrCurr), threshA);
%             end
%             kf(i, j, l).gammaMax = mean([kf_d(:,l).gammaMax]);
%             kf(i, j, l).sigmaAP = mean([kf_d(:,l).sigmaAP]);
%             kf(i, j, l).sigmaDP = mean([kf_d(:,l).sigmaDP]);
%         end

        for k = 1 : expNum
            pos = randi([1, lenEnv - lenFrame]);
            e = awgn(env(pos:pos+lenFrame), snrCurr, 'measured');
            kf_d(k) = KeyFeatures(e, threshA);
%             iteration = iteration + 1;
%             waitbar(iteration / cyclesNum);
        end
        kf(i, j).gammaMax = mean([kf_d.gammaMax]);
        kf(i, j).sigmaAP = mean([kf_d.sigmaAP]);
        kf(i, j).sigmaDP = mean([kf_d.sigmaDP]);
        kf(i, j).sigmaAF = mean([kf_d.sigmaAF]);
        kf(i, j).sigmaDF = mean([kf_d.sigmaDF]);
        kf(i, j).P = mean([kf_d.P]);
        ci(i, j).gammaMax = ConfInt([kf_d.gammaMax]', beta); 
        ci(i, j).sigmaAP  = ConfInt([kf_d.sigmaAP]', beta);
        ci(i, j).sigmaDP  = ConfInt([kf_d.sigmaDP]', beta);
        ci(i, j).sigmaAF  = ConfInt([kf_d.sigmaAF]', beta);
        ci(i, j).sigmaDF  = ConfInt([kf_d.sigmaDF]', beta);
        ci(i, j).P        = ConfInt([kf_d.P]', beta);
        ci(i, j).gammaMax = ci(i, j).gammaMax(:,1);
        ci(i, j).sigmaAP  = ci(i, j).sigmaAP(:,1);
        ci(i, j).sigmaDP  = ci(i, j).sigmaDP(:,1);
        ci(i, j).sigmaAF  = ci(i, j).sigmaAF(:,1);
        ci(i, j).sigmaDF  = ci(i, j).sigmaDF(:,1);
        ci(i, j).P        = ci(i, j).P(:,1);
    end
end
toc
% close(h);
disp('... done.');

%% Plot KFs
close all;
PlotKeyFeatures(kf, ci, snr, modNames);




