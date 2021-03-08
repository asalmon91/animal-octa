function [Gc, C1, C2] = lv_to_ml_dispCompFrame(in_frame, p, interpIndex, roi_xywh)
%user_assisted_dispComp asks the user to select an ROI for dispersion
%compensation

%% Imports
addpath(genpath('.'));

% % DEV/DB
% octa_ffname = 'F:\img\2019.11.27-DM_186302\OCTA\2019_11_27_OS\Raw\warmed\DM_186302-20191127_122057-OS.octa';
% scan = getScanObj(octa_ffname);
% in_frame = read_octa_frames(octa_ffname, scan, round(scan.B/2), false);
% % END DEV/DB

%% Optional inputs
if exist('p','var') == 0 || isempty(p) || ...
        exist('interpIndex', 'var') == 0 || isempty(interpIndex)
    p = 1:size(in_frame, 1);
    interpIndex = loadSpecCal();
end
if exist('roi_xywh', 'var') == 0 || isempty(roi_xywh) || ...
        numel(roi_xywh) ~= 4
    roi_xywh = [1, 50, size(in_frame, 2)-1, size(in_frame,1)/4];
end

in_frame = single(in_frame) - mean(in_frame, 2);
in_frame = resampleOCU(in_frame, p, interpIndex);

%% Optimize dispersion
k0 = p(end)/2;
C_vec = dispComp_fminbnd(in_frame, [], [], roi_xywh);
Gc = exp(1i*(C_vec(1)*(p-k0).^2 + C_vec(2)*(p-k0).^3));
C1 = C_vec(1);
C2 = C_vec(2);

% % DEV/DB
% fid = fopen('C:\Users\octa-user\Documents\MATLAB\dispComp.txt', 'w');
% fprintf(fid, 'C1: %1.6fe-05; C2: %0.6fe-09', C1*1e5, C2*1e9);
% fclose(fid);
% fft_frame = abs(fft(in_frame, [], 1));
% out_frame = abs(fft(in_frame .* Gc', [], 1));
% % out_frame = out_frame(1:size(in_frame,1)/2, :);
% figure;
% subplot(1,2,1);
% imagesc(fft_frame);
% title('Uncompensated');
% subplot(1,2,2);
% imagesc(out_frame);
% title('Compensated');
% % END DEV/DB

end

