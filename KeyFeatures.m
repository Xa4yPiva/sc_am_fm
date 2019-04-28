function kf = KeyFeatures(envelope, aThreshold)
%KEYFEATURES Summary of this function goes here
%   Detailed explanation goes here

len = length(envelope);
Nfft = 2 ^ nextpow2(len);

kf.P = CoeffP(envelope, Nfft);
[aNorm, aCN] = SubCoeffsA(envelope);
kf.gammaMax = CoeffGammaMax(aCN, Nfft);
[kf.sigmaAP, kf.sigmaDP] = CoeffSigma(envelope, aNorm, aThreshold);
[kf.sigmaAF, kf.sigmaDF] = CoeffSigmaF(envelope, aNorm, aThreshold);

end

