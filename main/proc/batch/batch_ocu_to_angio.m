%% Import
addpath(genpath('F:\OCT-software\animal-octa\main\proc'));
addpath('F:\OCT-software\animal-octa\development-testing\angio');

%% Get dataset
octa_path = uigetdir('F:\img', 'Select directory');
% octa_path = 'F:\img\2019.10.29-DM_180402\2019_10_29_OS\Raw';
octa_dir = dir(fullfile(octa_path, '*.tiff'));
% bg_path = strrep(octa_path, 'Raw', 'Calibration');
% bg_fname = 'bg-20191022_150043-OS.tiff';


%% Set up progress bar
wb = waitbar(0, sprintf('Processing %s%s...', ...
    octa_dir(1).name, octa_dir(1).name));
wb.Children.Title.Interpreter = 'none';
waitbar(0, wb, sprintf('Processing %s...', octa_dir(1).name));

%% Load spectrometer calibration
p = 1:2048;
k0 = p(end)/2;
interpIndex = loadSpecCal();

%% Set up SSADA filters
GM = split_spectrum(1:2048, 4, 0.2, 0);
nGM = size(GM, 1);

%% Start going through files
for ii=1:numel(octa_dir)
    % Get name
    octa_fname = octa_dir(ii).name;
    
    % Determine xB...
    % todo: change this into a structure with all the parameters
    [xB, fail] = getXB(fullfile(octa_path, octa_fname));
    if fail
        continue;
    end
    if xB < 2
        warning('Not an angio scan, only outputting amplitude');
    end
    
    % Prepare output
    if xB > 1
        angio_out_fname = strrep(octa_fname, '.tiff', '-ang.tiff');
        ra_out_fname = strrep(octa_fname, '.tiff', '-ra.tiff');
    end
    if xB > 2
        vista_out_fname = strrep(octa_fname, '.tiff', '-vista.tiff');
    end
    fft_out_name = strrep(octa_fname, '.tiff', '-fft.tiff');
    % todo: check that this file renaming worked or the original scans may
    % be overwritten
    
    % Get indexing info
    octa_info = imfinfo(fullfile(octa_path, octa_fname));
    fiv = 1:numel(octa_info); % Frame indexing vector
    fiv = reshape(fiv, xB, numel(octa_info)/xB);
    nB = size(fiv, 2);
    ht = octa_info(1).Height;
    wd = octa_info(2).Width;
    
    %% Optimize dispersion based on the middle frame
    % TODO: include frame selection then ROI selection
    % TODO: ask whether they want to calculate dispersion for just the
    % first scan and apply those parameters to the rest
    if exist('Gc', 'var') == 0
        mid_frame_index = round(numel(octa_info)/2);
        frame = single(imread(fullfile(octa_path, octa_dir(ii).name), ...
            mid_frame_index));
        frame = subtractBackground(frame);
        frame = resampleOCU(frame, p, interpIndex, wb);
        fft_frame = ocu_fft(frame, wb);

        % Get user-defined ROI
        f = figure;
        ax = gca;
        imagesc(fft_frame)
        title('Double click roi when done');
        dispCompROI = imrect(ax, ...
            [size(frame,2)/3, size(frame,1)/2/3, ...
            size(frame,2)/3, size(frame,1)/2/3]);
        roi = round(wait(dispCompROI));
        close(f);
        
        %% Optimize dispersion
        C_vec = dispComp_fminbnd(frame, [], [], roi);
        Gc = exp(1i*(C_vec(1)*(p-k0).^2 + C_vec(2)*(p-k0).^3));
        
    end
    
    %% Get background vector from whole volume
    bg = getBG(fullfile(octa_path, octa_dir(ii).name), wb);
    
    %% Read the scans as clusters
    for jj=1:nB
        frames = zeros(ht, wd, size(fiv, 1), 'single');
        for tt=1:xB
            frames(:,:,tt) = single(...
                imread(fullfile(octa_path, octa_fname), ...
                fiv(tt, jj)));
        end
        
        % Process
        frames = subtractBackground(frames, bg);
        frames = resampleOCU(frames, p, interpIndex);
        frames = applyDispComp(frames, Gc);
        
        % Write structural volume
        log_frames = ocu_fft(frames);
        log_frames = uint16(log10(log_frames+1).*25e3);
        for tt=1:xB
            if jj==1 && tt == 1
                wm = 'overwrite';
            else
                wm = 'append';
            end
            imwrite(log_frames(:,:,tt), ...
                fullfile(octa_path, fft_out_name), 'writemode', wm);
        end
        
        if xB > 1
            % Filter the spectrum
            filt_amp_frames = zeros(ht/2, wd, xB, nGM, 'single');
            for tt=1:xB
                for mm=1:nGM
                    filt_frame = frames(:,:,tt) .* GM(mm ,:)';
                    fft_frame = abs(fft(filt_frame, [] , 1))./ht;
                    filt_amp_frames(:,:,tt,mm) = fft_frame(1:ht/2, :);
                end
            end

            %% Register the unfiltered amp frames
            % Get unfiltered amp frames
            crop_roi = {25:ht/4, 15:wd}; % todo: figure out a better way
            reg_frames = zeros(ht/2, wd, xB, 'single');
            reg_roi_frames = zeros(...
                numel(crop_roi{1}), ...
                numel(crop_roi{2}), xB, 'single');
            for tt=1:xB
                fft_frame = abs(fft(frames(:,:,tt), [], 1))./ht;
                reg_frames(:,:,tt) = fft_frame(1:ht/2, :);
                reg_roi_frames(:,:,tt) = fft_frame(...
                    crop_roi{1}, crop_roi{2});
            end
            % Choose reference frame
            % for now, this is just the middle frame
            rfi = round(xB/2); % reference frame index
            ref_frame = reg_roi_frames(:,:,rfi);

            % Register amp frames
            reg_filt_amp_frames = filt_amp_frames;
            reg_full_frames = reg_frames;
            fprintf('Frame cluster %i\n', jj);
            for tt=1:xB
                if tt==rfi
                    % Don't register the reference frame to itself
                    continue;
                end
                % Register frame based on cropped unfiltered amplitude
                fixedRefObj = imref2d(size(ref_frame));
                mov_frame = reg_roi_frames(:,:,tt);
                movingRefObj = imref2d(size(mov_frame));

                % Phase correlation
                tform = imregcorr(mov_frame, movingRefObj, ...
                    ref_frame, fixedRefObj, ...
                    'transformtype', 'translation', 'Window', true);
                % Display results
                fprintf('dx: %0.2f, dy: %0.2f\n', ...
                    tform.T(3,1), tform.T(3,1));
                
                %% Register unfiltered full frames
                fixedRefObj  = imref2d(size(reg_full_frames(:,:,1)));
                movingRefObj = fixedRefObj;
                reg_full_frames(:,:,tt) = imwarp(...
                    reg_frames(:,:,tt), movingRefObj, ...
                    tform, 'OutputView', fixedRefObj, ...
                    'SmoothEdges', true);
                
                %% Register filtered frames
                % Reset spatial reference objects to filtered amp frames
                fixedRefObj  = imref2d(size(filt_amp_frames(:,:,1,1)));
                movingRefObj = fixedRefObj;
                for mm=1:nGM
                    % Apply transformation to filtered amp frames
                    reg_filt_amp_frames(:,:,tt,mm) = imwarp(...
                        filt_amp_frames(:,:,tt,mm), movingRefObj, ...
                        tform, 'OutputView', fixedRefObj, ...
                        'SmoothEdges', true);
                end
            end
            if jj==1
                wm = 'overwrite';
            else
                wm = 'append';
            end
            ra_frame = mean(reg_full_frames, 3);
            log_ra_frame = uint16(log10(ra_frame+1).*25e3);
            imwrite(log_ra_frame, ...
                fullfile(octa_path, ra_out_fname), 'writemode', wm);

            % Measure decorrelation on registered and filtered amp frames
            d = zeros(ht/2, wd, 'single');
            for tt=1:xB-1
                for mm=1:nGM
                    d = d + decorrelation(...
                        reg_filt_amp_frames(:,:,tt,mm), ...
                        reg_filt_amp_frames(:,:,tt+1,mm));
                end
            end
            octa_frame = 1 - d./nGM./(xB-1);

            % VISTA x2 (double TR)
            if xB > 2
                xTR = 2; % repetition time multiplier
                nComps = xB - xTR; % # comparisons
                vix = zeros(nComps, 2); % VISTA index
                for xx=1:nComps
                    vix(xx, :) = [xx, xx+xTR];
                end

                d = zeros(ht/2, wd, 'single');
                for tt=1:nComps
                    for mm=1:nGM
                        d = d + decorrelation(...
                            reg_filt_amp_frames(:,:,vix(tt,1),mm), ...
                            reg_filt_amp_frames(:,:,vix(tt,2),mm));
                    end
                end
                octa_frame_tr2 = 1 - d./nGM./(nComps);
                
            end

            % Threshold based on amplitude
            % todo: do this based on the whole volume
            avg_amp_frame = mean(reg_frames, 3);
            amp_mask = avg_amp_frame >= ...
                mean(avg_amp_frame(:)) + 2*std(avg_amp_frame(:));
%             d_mask = octa_frame >= ...
%                 mean(octa_frame(:)) + 0*std(octa_frame(:));
            octa_frame = octa_frame.*amp_mask;

            % Intensity transformation for writing to tiff
            octa_frame = uint16(octa_frame.*65535);

            if jj==1
                wm = 'overwrite';
            else
                wm = 'append';
            end
            imwrite(octa_frame, fullfile(octa_path, angio_out_fname), ...
                'writemode', wm);
            if xB > 2
                octa_frame_tr2 = uint16(octa_frame_tr2.*amp_mask.*65535);
                imwrite(octa_frame_tr2, fullfile(octa_path, vista_out_fname), ...
                    'writemode', wm);
            end
        end
        
        waitbar(jj/nB, wb);
    end
    
    waitbar(ii/numel(octa_dir), wb, ...
        sprintf('Done processing %s.', octa_fname));
end
close(wb);






