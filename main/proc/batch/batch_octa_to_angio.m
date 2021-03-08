%% Import
addpath(genpath('F:\OCT-software\animal-octa\main\proc'));
addpath('F:\OCT-software\animal-octa\main\fx\output\mat');
addpath('F:\OCT-software\animal-octa\development-testing\angio');

%% Get dataset
octa_path = uigetdir('F:\img', 'Select directory');
% octa_path = 'F:\img\2019.11.26-DM_175003\OCTA\2019_11_26_OS\Raw';
octa_dir = dir(fullfile(octa_path, '*.octa'));
out_path = strrep(octa_path, 'Raw', 'Processed');
if exist(out_path, 'dir') == 0
    mkdir(out_path);
end

%% Set up progress bar
wb = waitbar(0, sprintf('Processing %s%s...', ...
    octa_dir(1).name, octa_dir(1).name));
wb.Children.Title.Interpreter = 'none';
waitbar(0, wb, sprintf('Processing %s...', octa_dir(1).name));

%% Load spectrometer calibration
z = 2048;
p = 1:z;
k0 = p(end)/2;
interpIndex = loadSpecCal();

%% Set up SSADA filters
GM = split_spectrum(1:2048, 8, 0.2, 0);
nGM = size(GM, 1);

%% Set up parallel pool
current_pool = gcp('nocreate');
if isempty(current_pool)
    current_pool = parpool(8);
end

%% Start going through files
for ii=1:numel(octa_dir)
    % Get naming shortcuts
    octa_fname  = octa_dir(ii).name;
    octa_ffname = fullfile(octa_path, octa_fname);
    [~, octa_name, octa_ext] = fileparts(octa_ffname);
    
    % Determine xB...
    % todo: change this into a structure with all the parameters
    [scan, fail, err] = getScanObj(fullfile(octa_path, octa_fname));
    if fail
        disp(err);
        continue;
    end
    xB = scan.xB;
    xC = scan.C;
    if xB < 2 && xC < 2
        warning('Not an angio scan, only outputting amplitude');
    end
    
    % Prepare output
    angio_out_fname = [octa_name, '-ssada1.tiff'];
%     skip_ang = ...
%         exist(fullfile(out_path, angio_out_fname), 'file') ~= 0 | xB < 2;
    skip_ang = true;
    ra_out_fname = [octa_name, '-ra.tiff'];
    skip_ra = ...
        exist(fullfile(out_path, ra_out_fname), 'file') ~= 0 | ...
        (xB < 2 & xC < 2); 
    var_out_fname = [octa_name, '-var.tiff'];
    skip_var = exist(fullfile(out_path, var_out_fname), 'file') ~= 0 ...
        | (xB < 2 & xC < 2);
    vista_out_fname = [octa_name, '-ssada2.tiff'];
%     skip_vista = ...
%         exist(fullfile(out_path, vista_out_fname), 'file') ~= 0 | xB < 3;
    skip_vista = true;
    fft_out_fname = [octa_name, '-oct.tiff'];
%     skip_fft = exist(fullfile(out_path, fft_out_fname), 'file') ~= 0;
    skip_fft = true;
    reg_out_fname = [octa_name, '-reg.tiff'];
    skip_reg = true;
    skip_vec = [skip_ang, skip_ra, skip_reg, skip_var, skip_vista, skip_fft];
    if all(skip_vec)
        continue;
    end
    
    % Get indexing info
    B = scan.B;
    fiv = 1:B*xB*xC; % Frame indexing vector
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
    
    %% Optimize dispersion based on the middle frame
    % TODO: include frame selection then ROI selection
    % TODO: ask whether they want to calculate dispersion for just the
    % first scan and apply those parameters to the rest
    if exist('Gc', 'var') == 0
        mid_frame_index = round(B*xB/2);
        frame = single(read_octa_frames(octa_ffname, scan, mid_frame_index));
        frame = subtractBackground(frame);
        frame = resampleOCU(frame, p, interpIndex);
        fft_frame = ocu_fft(frame);

        % Get user-defined ROI
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
        C_vec = dispComp_fminbnd(frame, [], [], roi);
        Gc = exp(1i*(C_vec(1)*(p-k0).^2 + C_vec(2)*(p-k0).^3));
        
    end
    
    %% Get background vector from whole volume
%     tic
%     if exist('bg', 'var') == 0 || isempty(bg)
    bg = getBG(octa_ffname, scan, wb);
%     end
%     t = toc/B*xB;
    
    %% Read and process the scans as clusters
    if ~skip_ra
        ra_vol = zeros(ht/2, wd, B, 'uint16');
    end
    if ~skip_var
        var_vol = zeros(ht/2, wd, B, 'uint16');
    end
    
    
    parfor jj=1:B
        frame_indices = fiv(:, jj);
        frames = single(read_octa_frames(octa_ffname, scan, frame_indices));
        frames = subtractBackground(frames, bg);
        frames = resampleOCU(frames, p, interpIndex);
        frames = applyDispComp(frames, Gc);
        % mangitude image for registration
        fft_frames = abs(fft(frames,[],1));
        fft_frames = fft_frames(1:ht/2, :, :);
        
        % Write structural volume
        if ~skip_fft
            log_frames = uint16(log10(fft_frames./ht+1).*25e3);
            for tt=1:xRpt
                if jj==1 && tt == 1
                    wm = 'overwrite';
                else
                    wm = 'append';
                end
                imwrite(log_frames(:,:,tt), ...
                    fullfile(out_path, fft_out_fname), 'writemode', wm);
            end
        end
        
        if xRpt > 1
            % Filter the spectrum
            if ~skip_ang && ~skip_vista
                filt_amp_frames = zeros(ht/2, wd, xB, nGM, 'single');
                for tt=1:xRpt
                    for mm=1:nGM
                        filt_frame = frames(:,:,tt) .* GM(mm ,:)';
                        fft_frame = abs(fft(filt_frame, [] , 1))./ht;
                        filt_amp_frames(:,:,tt,mm) = fft_frame(1:ht/2, :);
                    end
                end
            end

            %% Register the unfiltered amp frames
            % Crop out background for registration
            crop_roi = {30:ht/4, 25:wd}; % todo: figure out a better way
            reg_frames = fft_frames;
            reg_roi_frames = reg_frames(crop_roi{1}, crop_roi{2}, :);
            % Choose reference frame
            % for now, this is just the middle frame
            rfi = round(xRpt/2); % reference frame index
            ref_frame = reg_roi_frames(:,:,rfi);

            % Register amp frames
            fprintf('Frame cluster %i\n', jj);
            for tt=1:xRpt
                if tt==rfi
                    % Don't register the reference frame to itself
                    continue;
                end
                % Register frame based on cropped unfiltered amplitude
                fixedRefObj = imref2d(size(ref_frame));
                mov_frame = reg_roi_frames(:,:,tt);
                movingRefObj = imref2d(size(mov_frame));

                % Phase correlation
                warn_msg = ''; %#ok<NASGU>
                tform = imregcorr(mov_frame, movingRefObj, ...
                    ref_frame, fixedRefObj, ...
                    'transformtype', 'translation', 'Window', true);
                warn_msg = lastwarn;
                dxdy = [tform.T(3,1), tform.T(3,1)];
                % Display results
                fprintf('dx: %0.2f, dy: %0.2f\n', dxdy(1), dxdy(2));
                if ~isempty(warn_msg)
                    dxdy = [0,0]; % Don't register if fit is bad
                end
                % Don't warp anything less than half a pixel
                if all(abs(dxdy) < 0.5)
                    continue;
                end
                
                %% Register unfiltered full frames
                fixedRefObj  = imref2d(size(reg_frames(:,:,1)));
                movingRefObj = fixedRefObj;
                reg_frames(:,:,tt) = imwarp(...
                    reg_frames(:,:,tt), movingRefObj, ...
                    tform, 'OutputView', fixedRefObj, ...
                    'SmoothEdges', true);
                
                % Output registered frames to check performance

                if ~skip_reg
                    if jj==1
                        wm = 'overwrite';
                    else
                        wm = 'append';
                    end
                    log_reg_frame = uint16(log10(reg_frames(:,:,tt)./ht+1).*25e3);
                    imwrite(log_reg_frame, ...
                        fullfile(out_path, reg_out_fname), 'writemode', wm);
                end
                
                %% Register filtered frames
                % Reset spatial reference objects to filtered amp frames
                if ~skip_ang && ~skip_vista
                    fixedRefObj  = imref2d(size(filt_amp_frames(:,:,1,1)));
                    movingRefObj = fixedRefObj;
                    for mm=1:nGM
                        % Apply transformation to filtered amp frames
                        filt_amp_frames(:,:,tt,mm) = imwarp(...
                            filt_amp_frames(:,:,tt,mm), movingRefObj, ...
                            tform, 'OutputView', fixedRefObj, ...
                            'SmoothEdges', true);
                    end
                end
            end
            
            if jj==1
                wm = 'overwrite';
            else
                wm = 'append';
            end
            if ~skip_ra
                ra_frame = mean(reg_frames, 3);
                log_ra_frame = uint16(log10(ra_frame./ht+1).*25e3);
                ra_vol(:,:,jj) = log_ra_frame;
%                 imwrite(log_ra_frame, ...
%                     fullfile(out_path, ra_out_fname), 'writemode', wm);
            end
            if ~skip_var
                var_frame = var(reg_frames, [], 3) ./ mean(reg_frames,3);
                var_frame = var_frame./max(var_frame(:));
                c_range = [0, mean(var_frame(:)) + 3*std(var_frame(:))];
                var_frame = imadjust(var_frame, c_range);
                var_frame = uint16(var_frame.*(2^16-1));
                var_vol(:,:,jj) = var_frame;
%                 imwrite(var_frame, ...
%                     fullfile(out_path, var_out_fname), 'writemode', wm);
            end
            
            % Measure decorrelation on registered and filtered amp frames
            if ~skip_ang
                d = zeros(ht/2, wd, 'single');
                for tt=1:xRpt-1
                    for mm=1:nGM
                        d = d + decorrelation(...
                            filt_amp_frames(:,:,tt,mm), ...
                            filt_amp_frames(:,:,tt+1,mm));
                    end
                end
                octa_frame = 1 - d./nGM./(xRpt-1);
            end

            % VISTA x2 (double TR)
            if xRpt > 2 && ~skip_vista
                xTR = 2; % repetition time multiplier
                nComps = xRpt - xTR; % # comparisons
                vix = zeros(nComps, 2); % VISTA index
                for xx=1:nComps
                    vix(xx, :) = [xx, xx+xTR];
                end

                d = zeros(ht/2, wd, 'single');
                for tt=1:nComps
                    for mm=1:nGM
                        d = d + decorrelation(...
                            filt_amp_frames(:,:,vix(tt,1),mm), ...
                            filt_amp_frames(:,:,vix(tt,2),mm));
                    end
                end
                octa_frame_tr2 = 1 - d./nGM./(nComps);
            end

            % Threshold based on amplitude
            % todo: do this based on the whole volume
            if ~skip_ang && ~skip_vista
                avg_amp_frame = mean(reg_frames, 3);
                amp_mask = avg_amp_frame >= ...
                    mean(avg_amp_frame(:)) + std(avg_amp_frame(:));
    %             d_mask = octa_frame >= ...
    %                 mean(octa_frame(:)) + 0*std(octa_frame(:));
                octa_frame = octa_frame.*amp_mask;

                % Intensity transformation for writing to tiff
                octa_frame = uint16(octa_frame.*(2^16-1));
            end

            if jj==1
                wm = 'overwrite';
            else
                wm = 'append';
            end
            if ~skip_ang
                imwrite(octa_frame, fullfile(out_path, angio_out_fname), ...
                    'writemode', wm);
            end
            if xRpt > 2 && ~skip_vista
                octa_frame_tr2 = uint16(octa_frame_tr2.*amp_mask.*(2^16-1));
                imwrite(octa_frame_tr2, fullfile(out_path, vista_out_fname), ...
                    'writemode', wm);
            end
        end
        
%         waitbar(jj/B, wb, 'Processing frames');
    end
    
    for jj=1:B
        if jj==1
            wm = 'overwrite';
        else
            wm = 'append';
        end
        if ~skip_ra
            imwrite(ra_vol(:,:,jj), fullfile(out_path, ra_out_fname), ...
                'writemode', wm);
        end
        if ~skip_var
            imwrite(var_vol(:,:,jj), fullfile(out_path, var_out_fname), ...
                'writemode', wm);
        end
        
        waitbar(jj/B, wb, 'Writing videos');
    end
    
    waitbar(ii/numel(octa_dir), wb, ...
        sprintf('Done processing %s.', octa_fname));
end
close(wb);






