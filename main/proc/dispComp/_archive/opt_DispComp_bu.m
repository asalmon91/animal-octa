function [inv_sharp, IdFa] = opt_DispComp(CX, ...
    img, new_index, IdFa, xm, sharps, isC2, CY)
%opt_DispComp optimizes the dispersion compensation coefficients
%   todo: write detailed description

%% Get simple indexing variables
n=size(img,1);
p=1:n;
k0=n/2;

%% Get dispersion compensation vector
if isC2
    Gc = exp(1i*(CX*(p-k0).^2 + CY*(p-k0).^3));
else
    Gc = exp(1i*(CY*(p-k0).^2 + CX*(p-k0).^3));
end

%% Sharpness
% sharps = zeros(size(img,2), 1);
for ii=1:size(img,2)
    % Flat field correction for reduction of fixed noise pattern
    x1 = (((double(img(:,ii)))-xm));

    % Resampling data in linearized K-space with Spline interpolation
    IdiL = interp1(p, x1, new_index, 'spline');
    % Apply dispersion compensation
    IdiLc=IdiL.*Gc;
    
    % FFT
    IdF = abs(fft(IdiLc));
    IdFa(:,ii)=IdF(1:n/2); %Build the image

    sharps(ii) = getSharpness(IdFa(:,ii));
end

inv_sharp = 1/mean(sharps);




end

