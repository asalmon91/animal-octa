function [oct, Gc, svp, times] = proc_oct_gpu(octa_ffname, Gc, profiling)
%proc_oct_gpu Processes an OCT-A volume on the GPU

%% Imports
addpath(genpath('.'));

%% Reset GPU
gpuDevice(1);

%% Load spectrometer calibration
z = 2048;
p = 1:z;
interpIndex = loadSpecCal();

%% Get scan object
[scan, fail, err] = getScanObj(octa_ffname);
if fail
    error(err);
end

%% Get dispersion compensation
if exist('Gc', 'var') == 0 || isempty(Gc)
    Gc = user_assisted_dispComp(octa_ffname, scan, p, interpIndex);
end

%% Get indexing info
B = scan.B;
xB = scan.xB;
xC = scan.C;
N = B*xB*xC;
fiv = 1:N; % Frame indexing vector
if xB > 1
    fiv = reshape(fiv, xB, B);
    xRpt = xB;
elseif xC > 1
    fiv = reshape(fiv, B, xC)';
    xRpt = xC;
else
    xRpt = 1;
end
ht = z;
wd = scan.A;

%% Preallocate a gpuArray for each cluster of frames
raw_frames      = gpuArray(zeros(           ht,    wd,         xRpt,    'uint16'));
proc_frames     = gpuArray(complex(zeros(   ht,    wd*xRpt,    1,       'single')));
fft_frames      = gpuArray(zeros(           ht/2,  wd*xRpt,    1,       'single'));
log_fft_frames  = gpuArray(zeros(           ht/2,  wd*xRpt,    1,       'uint16'));
% And a CPU array for storing the final product
oct = zeros(ht/2, wd, N, 'uint16');

%% Profile
if exist('profiling', 'var') == 0 || isempty(profiling)
    profiling = false;
end
if profiling
    proc_modules = {'Read', 'BG', 'InterpDispCompFFT', 'Int', 'ToHost'};
    read_c      = strcmp(proc_modules, 'Read');
    bg_c        = strcmp(proc_modules, 'BG');
    proc_c      = strcmp(proc_modules, 'InterpDispCompFFT');
    int_c       = strcmp(proc_modules, 'Int');
    transfer_c  = strcmp(proc_modules, 'ToHost');
    times = zeros(B, numel(proc_modules));
end

%% Read and Process
dims = [ht, wd, N];
fid  = fopen(octa_ffname, 'r');
tic
for frame_cluster_i = 1:B
    % Get frame indices
    frames_in_this_cluster = fiv(:, frame_cluster_i);
    
    tic
    %% Frame Cluster Loop
    for sub_frame_i = 1:numel(frames_in_this_cluster)
        %% Read
        % Move to start position of this frame
        status = fseek(fid, ...
            (prod(dims(1:2))*(frames_in_this_cluster(sub_frame_i)-1))*2, ... 
            'bof'); % x2 16-bit = 2bytes
        if status ~= 0
            error('I don''t know binary as well as I thought');
        end

        % Read frame
        raw_frames(:, :, sub_frame_i) = ...
            fread(fid, [dims(1), dims(2)], 'uint16', 0, 'l');
    end
    if profiling
        times(frame_cluster_i, read_c) = toc;
    end
    
    %% DC Mitigation (background subtraction)
    tic
    bg = single(mean(raw_frames, 2));
    % Stack frames next to each other for use with interp1
    proc_frames(:,:) = reshape(raw_frames, [ht, wd*xRpt, 1]);
    proc_frames(:,:) = proc_frames(:,:) - bg;
    if profiling
        times(frame_cluster_i, bg_c) = toc;
    end
    
    %% Resampling, Dispersion Compensation, and FFT
    tic
    proc_frames(:,:) = fft(...
        interp1(...
        p', proc_frames, interpIndex', 'linear') .* ...
        Gc', [] , 1);
    fft_frames(:,:) = abs(proc_frames(1:ht/2, :))./ht;
    if profiling
        times(frame_cluster_i, proc_c) = toc;
    end
    
    %% Intensity transformation for human-viewing
    tic
    log_fft_frames(:,:) = uint16(log10(fft_frames+1).*25e3);
    if profiling
        times(frame_cluster_i, int_c) = toc;
    end
    
    %% Transfer to CPU as individual frames
    tic
    oct(:,:, frames_in_this_cluster) = ...
        gather(reshape(log_fft_frames, [ht/2, wd, xRpt]));
    if profiling
        times(frame_cluster_i, transfer_c) = toc;
    end
    
    % Simple progress update
    if mod(frame_cluster_i, 100) == 0
        fprintf('%0.2f%s', frame_cluster_i/B*100, '%');
        if frame_cluster_i==B
            fprintf('\n');
        else
            fprintf(', ');
        end
    end
end
fclose(fid);
t = toc;
fprintf('Total time: %0.2fs\n', t);
fprintf('Time/frame: %0.2fs\n', t/N);

%% Get SVP
svp = squeeze(mean(oct, 1))';
svp = svp-min(svp(:));
svp = uint8(svp./max(svp(:)).*255);

%% Output
svp_out_ffname = strrep(octa_ffname, '.octa', '.tif');
if strcmp(svp_out_ffname, octa_ffname)
    error('Failed to rename %s', octa_ffname);
end
imwrite(svp, svp_out_ffname);

if ~profiling
    times = [];
end
% % DEV/DB
% figure;
% imagesc(svp);
% figure;
% imshow(svp);
% % END DEV/DB

% % DEV/DB
% figure;
% for ii=1:N
%     imshow(oct(:,:,ii));
%     pause(1/30);
% end
% % END DEV/DB

% DEV/DB
figure;
hold on;
for ii=1:size(times, 2)
%     histogram(times(:, ii).*1e3);
    [y,x] = histcounts(times(:, ii).*1e3, N/10);
    stem(x(1:end-1),y);
end
legend(proc_modules, 'location','northeastoutside');
hold off;
grid on
set(gca,'xscale','log', 'yscale','log','tickdir','out');
xlabel('Time (ms)');
ylabel('# Frames');
ylim([1e-1, N])
% END DEV/DB


end

