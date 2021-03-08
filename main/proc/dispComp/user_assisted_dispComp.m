function Gc = user_assisted_dispComp(oct_ffname, scanObj, p, interpIndex, frame)
%user_assisted_dispComp asks the user to select an ROI for dispersion
%compensation

%% Get middle frame
% (most likely to be one of the better frames)
if exist('frame', 'var') == 0 || isempty(frame)
    mid_frame_index = round(scanObj.B*scanObj.xB/2);
    % Read and process
    frame = read_octa_frames(oct_ffname, scanObj, mid_frame_index);
else
    if isa(frame, 'gpuArray')
        frame = gather(frame);
    end
end
frame = subtractBackground(single(frame));
frame = resampleOCU(frame, p, interpIndex);
fft_frame = ocu_fft(frame);

%% Get user-defined ROI
f = figure;
ax = gca;
imagesc(fft_frame)
title('Double click roi when done');
dispCompROI = imrect(ax, ...
    [size(fft_frame,2)/3, size(fft_frame,1)/2/3, ...
    size(fft_frame,2)/3, size(fft_frame,1)/2/3]);
roi = round(wait(dispCompROI));
close(f);

%% Optimize dispersion
k0 = p(end)/2;
C_vec = dispComp_fminbnd(frame, [], [], roi);
Gc = exp(1i*(C_vec(1)*(p-k0).^2 + C_vec(2)*(p-k0).^3));

% % DEV/DB
% fft_frame = abs(fft(frame .* Gc', [], 1));
% fft_frame = fft_frame(1:1024, :);
% figure;
% imagesc(fft_frame);
% % END DEV/DB

end

