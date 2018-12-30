function P = CoeffP(envelope, Nfft)
%COEFFP Summary of this function goes here
%   Detailed explanation goes here

spectrumAbs = abs(fftshift(fft(envelope, Nfft))) / (Nfft/2);
f0 = Nfft / 2 + 1; % Check this
Pl = sum(spectrumAbs(1 : f0 - 1)   .* spectrumAbs(1 : f0 - 1));
Pu = sum(spectrumAbs(f0 + 1 : end) .* spectrumAbs(f0 + 1 : end));
P = (Pl - Pu) / (Pl + Pu);

end

