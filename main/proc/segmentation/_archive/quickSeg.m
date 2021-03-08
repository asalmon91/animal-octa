function [top, bottom] = quickSeg(vol)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



% vol = single(vol)-bg;
% vol = resampleOCU(vol, p, interpIndex);
% vol = applyDispComp(vol, Gc);
% vol = ocu_fft(vol);

vol_roi = vol(30:end, 15:end, :);
avg_vol_roi = zeros(size(vol_roi,1), size(vol_roi,2), size(vol_roi,3)/8, 'single');
for ii=1:size(avg_vol_roi, 3)
    avg_vol_roi(:,:,ii) = mean(vol_roi(:,:,(1:8)+(8*(ii-1))), 3);
end

w = gausswin(2048/10);
blur_vol = diff(filter(w,1,avg_vol_roi,[],1), 1, 1);
diff_vol = abs(diff(avg_vol_roi.^3,[],1));
[~, maxIndex] = max(blur_vol, [], 1);
[~, minIndex] = min(blur_vol, [], 1);
% Account for missing element due to diff
maxIndex = maxIndex - 1;
minIndex = minIndex - 1;

mfs = 7;
top = medfilt2(squeeze(maxIndex), [mfs,mfs]);
bottom = medfilt2(squeeze(minIndex), [mfs,mfs]);
top(top<0) = 0;
bottom(bottom<0) = 0;
% Could certainly be improved

figure;
svp = sum(avg_vol_roi(round(top+1:top+200), :, :), 1);
imagesc(squeeze(svp)');

figure;
for ii=1:size(avg_vol_roi, 3)
    imagesc(avg_vol_roi(:,:,ii));
    hold on;
    plot(1:size(avg_vol_roi, 2), top(ii,:));
    plot(1:size(avg_vol_roi, 2), bottom(ii,:));
    hold off;
    
    pause(0.2);
end

end

