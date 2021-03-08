function newIndex = loadSpecCal()
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Spectrometer calibration
SC0 = 7.52368E+02;
SC1 = 9.89288E-02;
SC2 = -5.15442E-06;
SC3 = -1.65648E-10;
p = 1:2048;
n = max(p);
lambda_nm = SC0 + SC1*p + SC2*p.^2 + SC3*p.^3;

%% Conversion to linear k-space
% The non-linear k-space
K = 2*pi./lambda_nm;
% The rearranged sampling interval in the linear K-space
DeltaK = abs(K(1) - K(end))/n;
% The linear K-space
Klin = flip(K(n) + p * DeltaK);
newIndex = interp1(K,p,Klin);

end

