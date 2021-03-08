%% Spectrometer Calibration Function
C0 = 7.52368E+02;
C1 = 9.89288E-02;
C2 = -5.15442E-06;
C3 = -1.65648E-10;
p = 1:2048;
n = max(p);
lambda_nm = C0 + C1*p + C2*p.^2 + C3*p.^3;

%% Conversion to linear k-space
% The non-linear k-space
K = 2*pi./lambda_nm;
% The rearranged sampling interval in the linear K-space
DeltaK = abs(K(1) - K(end))/n;
% The linear K-space
Klin = flip(K(n) + p * DeltaK);

newIndx = interp1(K,p,Klin);
