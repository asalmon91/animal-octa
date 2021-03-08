function new_Indx = LambdaToK_Standard4_as(I)

%% Created by Farid Atry
% Modified by Alex Salmon - 2018.10.16 - axexsalmon@gmail.com
if exist('I', 'var') == 0
    [img_fname, img_path] = uigetfile('*.tif', 'Select interference image');
    img = imread(fullfile(img_path, img_fname));
    I = mean(img, 1);
end

%% Constants
HP = 0.0075; % low pass
LP = 1-HP; % high pass
TAPER = 0.5;

%% clean the calibration data
% Bandpass filter to remove source spectrum;
I_FFT = fft(I,[],2);
I_FFT = abs(I_FFT(:,1:ceil(size(I_FFT,2)/2)));
I_FFT(:, 1:ceil(HP * size(I_FFT, 2))) = 0;
I_FFT(:, ceil(LP * size(I_FFT, 2)):end) = 0;
% Interference becomes the strongest peak in the Fourier domain signal
[~, Indx] = max(I_FFT,[],2);

% Filter signal with Tukey Window
Tukey_mask = zeros(size(I,1), size(I,2)/2+1);
for ii=1:size(I,1)
    x = ceil(Indx(ii)/2);
%     if Indx(ii) > x && (Indx(ii)+2*x) < size(I,2)/2
    Tukey_mask(ii, Indx(ii)-x:Indx(ii)+2*x-1) = tukeywin(3*x, TAPER)';
%     elseif Indx(ii) < x
%         Tukey_mask(ii,1:Indx(ii)+3*x-1) = tukeywin(3*x,.2)';
%     else
%         Tukey_mask(ii,end-350+1:end) = tukeywin(350,.2)';
%     end
end
% Apply filter
I_FFT = fft(I,[],2).*[Tukey_mask, Tukey_mask(:,end-1:-1:2)];



%% Hilbert transform to exract linear phase samples
A = fftshift(fft(ifft(I_FFT,[],2),[],2),2);
A(:, 1:fix(size(A, 2)/2)+1) = 0;
A = ifftshift(A, 2);
b = ifft(A,[],2);
B = angle(b);
B = unwrap(B,[],2);
B = mean(B,1);

% n = 5;
% p = polyfit(...
%     (ceil(HP * length(I)):ceil(LP * length(I))), ...
%     B(ceil(HP * length(I)):ceil(LP * length(I))), n);
% Indx = 1:length(I);
% ExtendedPhase2 = p(n+1);
% for ii = 1:size(p,2)-1
%     ExtendedPhase2 = ExtendedPhase2 + p(ii)*Indx.^(n-ii+1);
% end
%% Filter out ends, extend to original limits
f5 = fit(...
    (ceil(HP * length(I)):ceil(LP * length(I)))', ...
    B(ceil(HP * length(I)):ceil(LP * length(I)))', 'poly5');
ExtendedPhase = feval(f5, 1:size(I, 2));

%% Obtain equispaced, monotonically increasing phasor
% UniPhase = ExtendedPhase(1): ...
%     (ExtendedPhase(end)-ExtendedPhase(1))/(size(I, 2)-1): ... % step size
%     ExtendedPhase(end);
UniPhase2 = linspace(ExtendedPhase(1), ExtendedPhase(end), size(I, 2));
new_Indx = spline(ExtendedPhase, 1:size(I, 2), UniPhase2);

%% Display results
figure
subplot(3,1,1);
plot(mean(I,1), 'k');
axis tight
title('Averaged spectrometer signal');
xlabel('Camera pixel');
ylabel('Intensity');
set(gca, 'tickdir', 'out', 'box', 'off');

subplot(3,1,2)
title('Frequency analysis');
hold on
plot(abs(fft(I(1,:))), 'k') % Original
plot(abs(I_FFT(1,:)), '--r') % Filtered
hold off
xlabel('Frequency (AU)');
ylabel('Power (AU)');
legend({'Original', 'Filtered'}, 'location', 'northeast');
axis tight
set(gca, 'tickdir', 'out');

subplot(3,1,3)
title('Unwrapped phase');
hold on
plot(B, 'k') % Measured
plot(ExtendedPhase,'--r') % Fit
hold off
xlabel('Camera pixel');
ylabel('Phase (rad)');
legend({'Measured', 'Fit'}, 'location', 'southeast');
axis tight
set(gca, 'tickdir', 'out');



