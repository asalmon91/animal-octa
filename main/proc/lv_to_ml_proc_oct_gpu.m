function [mvp_out_ffname, times, Gc] = lv_to_ml_proc_oct_gpu(octa_ffname, Gc, quickFeedBack)
%proc_oct_gpu Processes an OCT-A volume on the GPU

%% Imports
addpath(genpath('.'));

% DEV/DB
% octa_ffname = 'F:\img\2019.11.27-DM_164807\2019_11_27_OS\Raw\warmed\DM_164807-20191127_203857-OS.octa';
% END DEV/DB

%% Constants
BIT8 = 2^8-1;
BIT16 = 2^16-1;
profiling = true;

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
raw_frames      = gpuArray(zeros(           ht,    wd,          xRpt,   'uint16'));
proc_frames     = gpuArray(complex(zeros(   ht,    wd*xRpt,     1,      'single')));
fft_frames      = gpuArray(zeros(           ht/2,  wd*xRpt,     1,      'single'));
reg_frames      = gpuArray(zeros(           ht/2,  wd,          xRpt,   'single'));
log_fft_frames  = gpuArray(zeros(           ht/2,  wd,          1,      'uint16'));
octa_frame      = gpuArray(zeros(           ht/2,  wd,          1,      'uint16'));
% And a CPU array for storing the final product
oct     = zeros(ht/2, wd, B, 'uint16'); % Log-transformed
octa    = zeros(ht/2, wd, B, 'uint16');

%% Profile
if exist('profiling', 'var') == 0 || isempty(profiling)
    profiling = false;
end
if profiling
    proc_modules = {'Read', 'BG', 'InterpDispCompFFT', 'Reg', 'Angio', 'ToHost'};
    read_c      = strcmp(proc_modules, 'Read');
    bg_c        = strcmp(proc_modules, 'BG');
    proc_c      = strcmp(proc_modules, 'InterpDispCompFFT');
    reg_c       = strcmp(proc_modules, 'Reg');
    angio_c     = strcmp(proc_modules, 'Angio');
    transfer_c  = strcmp(proc_modules, 'ToHost');
    times = zeros(B, numel(proc_modules));
end

%% Read and Process
dims = [ht, wd, N];
try
    for frame_cluster_i = 1:B
        % Get frame indices
        frames_in_this_cluster = fiv(:, frame_cluster_i);
        
        tic
        %% Frame Cluster Loop
        fid  = fopen(octa_ffname, 'r');
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
        fclose(fid);
        if profiling
            times(frame_cluster_i, read_c) = toc;
        end
        
        %% DC Mitigation (background subtraction)
        tic
        %     bg = single(mean(raw_frames, 2));
        % Stack frames next to each other for use with interp1
        proc_frames(:,:) = reshape(raw_frames, [ht, wd*xRpt, 1]);
        bg = mean(proc_frames, 2);
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
        
        if xRpt > 1
            %% Register frames
            start_reg_time = datetime('now'); % tic toc doesn't work here
            reg_frames(:,:,:) = reshape(fft_frames, [ht/2, wd, xRpt]);
            if ~quickFeedBack
                reg_frames(:,:,:) = regFrames(reg_frames, 'accurate');
            end
            
            %% Registered-average
            if quickFeedBack
                oct_frame = reg_frames(:,:,1);
            else
                oct_frame = mean(reg_frames, 3);
            end
            if profiling
                times(frame_cluster_i, reg_c) = ...
                    seconds(datetime('now') - start_reg_time);
            end
            
            %% Get angiogram
            tic
            octa_frame(:,:) = uint16(get_fsada(reg_frames).*BIT16);
            if profiling
                times(frame_cluster_i, angio_c) = toc;
            end
        else
            oct_frame = fft_frames;
        end
        
        %% Intensity transformation for human-viewing
        tic
        log_fft_frames(:,:) = uint16(log10(oct_frame+1).*25e3);
        if profiling
            times(frame_cluster_i, reg_c) = toc;
        end
        
        %% Transfer to CPU as individual frames
        tic
        % Structural
        if isa(log_fft_frames, 'gpuArray')
            oct(:,:, frame_cluster_i) = gather(log_fft_frames);
        else
            oct(:,:, frame_cluster_i) = log_fft_frames;
        end
        % Angio
        if xRpt > 1
            if isa(octa_frame, 'gpuArray')
                octa(:,:, frame_cluster_i) = gather(octa_frame);
            else
                octa(:,:, frame_cluster_i) = octa_frame;
            end
        end
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
catch
    warning('Something failed');
    fclose(fid);
end

%% Get MVP
mvp = squeeze(mean(oct, 1))';
mvp = mvp-min(mvp(:));
mvp = uint8(mvp./max(mvp(:)).*255);

%% Output MVP
mvp_out_ffname = strrep(octa_ffname, '.octa', '.tif');
if strcmp(mvp_out_ffname, octa_ffname)
    error('Failed to rename %s', octa_ffname);
end
imwrite(mvp, mvp_out_ffname);

%% Quickly segment scan and get an angiogram
if xRpt > 1
%     raw_surface_seg = zeros(size(octa,3), size(octa,2));
%     for ii=1:B
%         raw_surface_seg(ii,:) = quickSegSurface(single(octa(:,:,ii)));
%     end
%     surface_seg = slowScanSmooth(raw_surface_seg);
%     
%     scp_yh = [10, 15];
%     mcp_yh = [25, 20];
%     dcp_yh = [45, 20];
%     cap_yh = [scp_yh; mcp_yh; dcp_yh];
%     angiograms = zeros([size(surface_seg), 3]);
%     for ii=1:3
%         angiograms(:,:,ii) = ...
%             squeeze(sum(octa(round(...
%             surface_seg+cap_yh(ii,1) : ...
%             surface_seg+cap_yh(ii,1)+cap_yh(ii,2)), :, :), 1))';
%     end
    figure; imagesc(squeeze(mean(octa, 1))');
    
%     % DEV/DB
%     if quickFeedBack
%         figure;
%         for ii=1:size(angiograms, 3)
%             subplot(1,size(angiograms, 3),ii)
%             imagesc(angiograms(:,:,ii));
%             axis square
%         end
%     end
%     % % END DEV/DB
end

%% Write out oct_vv-compatible .avi's
% Gather file information
% if ~quickFeedBack
[in_path, in_name, ~] = fileparts(octa_ffname);

out_path = fullfile(in_path, '..', 'Processed');
if exist(out_path, 'dir') == 0
    mkdir(out_path);
end

% Get output file names
oct_out_ffname  = fullfile(out_path, [in_name, '-ra.avi']);
octa_out_ffname = fullfile(out_path, [in_name, '-ada1.avi']);

% Contrast stretch
cs_range = stretchlim(oct(:));
for ii=1:size(oct, 3)
    oct(:,:,ii) = imadjust(oct(:,:,ii), cs_range);
end

% Move each to 8-bit
oct_vol_out  = uint8(single(oct)./BIT16*BIT8);
octa_vol_out = uint8(single(octa)./BIT16*BIT8);

% Crop range
% octa_vol_out = octa_vol_out(crop(1):crop(2), :, :); % crop DC + Autocorrelation
% octa_vol_out = flip(octa_vol_out, 3); % flip Z
OCX_2_AVI(oct_vol_out,  oct_out_ffname);
if xRpt > 1
    OCX_2_AVI(octa_vol_out, octa_out_ffname);
end
% end

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

% % DEV/DB
% figure;
% hold on;
% for ii=1:size(times, 2)
% %     histogram(times(:, ii).*1e3);
%     [y,x] = histcounts(times(:, ii).*1e3, round(N/10));
%     stem(x(1:end-1),y);
% end
% legend(proc_modules, 'location','northeastoutside');
% hold off;
% grid on
% set(gca,'xscale','log', 'yscale','log','tickdir','out');
% xlabel('Time (ms)');
% ylabel('# Frames');
% ylim([1e-1, N])
% % END DEV/DB

% % % DEV/DB
% figure;
% subplot(1,2,1);
% surf(raw_surface_seg), shading flat;
% subplot(1,2,2);
% surf(surface_seg), shading flat;
% % END DEV/DB

end

