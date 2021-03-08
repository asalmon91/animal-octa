function ocu_mat = ocu_gpu_fft(ocu_mat)
%ocu_fft Calculates the FFT of the resampled, compensated spectrometer
%image

% Get dimensions
[nZ, ~, ~] = size(ocu_mat);

% Do the fft
ocu_mat = abs(fft(ocu_mat, [], 1));
ocu_mat = ocu_mat(1:nZ/2, :, :)./nZ;

end

