%% Reset GPU
gpuDevice(1);

%% Load an OCT-A scan
octa_ffname = 'F:\img\2019.11.26-DM_175003\OCTA\2019_11_26_OS\Raw\DM_175003-20191126_153748-OS.octa';
scan = getScanObj(octa_ffname);

%% Get frame indices
N = scan.B*scan.xB;
fiv = reshape(1:N, scan.xB, scan.B);

%% Calibration vectors
ht = 2048;
p = 1:ht;
interpIndex = loadSpecCal();
Gc = [];

%% Preallocate storage for output
oct_vol     = zeros(1024, scan.A, scan.B, 'single');
octa_vol    = uint16(oct_vol);

%% This might take a while
wb = waitbar(0, sprintf('Reading %s...', octa_ffname));

%% Profiling
times = nan(scan.B, 1);
for ii=1:scan.B
    frame_idx = fiv(:,ii);
    in_frames = read_octa_frames(octa_ffname, scan, frame_idx, false);
    % Stack side-by-side because labview can only accept 2D arrays
    in_frames = reshape(in_frames, [size(in_frames, 1), scan.A*scan.xB, 1]);
    
    %% Run processing
    tic
    [oct_vol(:,:,ii), octa_vol(:,:,ii), Gc] = proc_buffered_oct_gpu(...
        in_frames, scan.A, Gc, p, interpIndex);
    times(ii) = toc;
    
    waitbar(ii/scan.B, wb, sprintf('Mean proc time: %0.3fs', ...
        mean(times, 'omitnan')))
end
close(wb);

%% Get SVPs
% Structural
svp = squeeze(mean(oct_vol, 1))';
svp = svp-min(svp(:));
svp = uint8(svp./max(svp(:)).*255);

% Angio
octa = squeeze(mean(octa_vol, 1))';
octa = octa-min(octa(:));
octa = uint8(octa./max(octa(:)).*255);

%% Display results
figure;
% Profile
subplot(2,2,1);
histogram(times, round(scan.B/10))
set(gca, 'xscale', 'log', 'yscale', 'log', 'ylim', [1e-1, scan.B])
xlabel('Time (s)');
ylabel('# Frames');
axis square
title(sprintf('Mean proc time: %0.3fs/frame', mean(times)));

% Structural SVP
subplot(2,2,2);
imshow(svp);

% Angio SVP
subplot(2,2,3);
imshow(octa);

