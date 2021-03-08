function [realImg, complexImg] = dispComp(camImg, dispCoeff1, dispCoeff2)
%dispComp applies the dispersion compensation formula to camImg
%   at this point camImg has already had the background subtracted and has
%   been resampled by the linear k-space vector.
%   camImg is a m*n U16 matrix where m == #A-scans/B-scan and n is the
%   number of pixels in the camera

%% Simple indexing variables
n=size(camImg,2);
p=1:n;
k0=n/2;

%% Get dispersion compensation vector
Gc = exp(1i*(dispCoeff1*(p-k0).^2 + dispCoeff2*(p-k0).^3));

%% Dispersion compensation and FFT
realImg = zeros(size(camImg,1), k0);
complexImg = complex(realImg);
for ii=1:size(camImg, 1)
    % FFT
    fft_camImg          = fft(camImg(ii,:).*Gc);
    complexImg(ii,:)    = fft_camImg(1:k0);
    realImg             = abs(complexImg);
    realImg(ii,:)       = realImg(1:k0);
end

end

