function oct = ocu_fft(ocu_mat, wb)
%ocu_fft Calculates the FFT of the resampled, compensated spectrometer
%image

[nZ, nA, nFrames] = size(ocu_mat);
k0 = nZ/2;

oct = zeros(k0, nA, nFrames, 'single');
for ii=1:size(ocu_mat, 3)
    mag_oct = abs(fft(ocu_mat(:,:,ii),[],1));
    oct(:,:,ii) = single(mag_oct(1:k0, :));
    
    if mod(ii, 10) == 0 && exist('wb', 'var') ~= 0 && ~isempty(wb)
        waitbar(ii/size(ocu_mat, 3), wb, 'Calculating FFT');
    end
end




end

