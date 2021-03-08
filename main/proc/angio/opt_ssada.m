%% Import
addpath(genpath('F:\OCT-software\animal-octa\main\proc'));
addpath('F:\OCT-software\animal-octa\development-testing\angio');

%% Analysis Constants
NSPLIT  = 1:15;
NBW     = 0.1:0.2:2.0;

%% Get dataset
[octa_fname, octa_path] = uigetfile('F:\img', 'Select directory', ...
    'multiselect', 'off');

%% Set up progress bar
wb = waitbar(0, sprintf('Processing %s%s...', ...
    octa_fname, octa_fname));
wb.Children.Title.Interpreter = 'none';
waitbar(0, wb, sprintf('Processing %s...', octa_fname));

%% Load spectrometer calibration
p = 1:2048;
k0 = p(end)/2;
interpIndex = loadSpecCal();

%% Determine xB...
% todo: change this into a structure with all the parameters
[xB, fail] = getXB(fullfile(octa_path, octa_fname));
if fail || xB < 2
    return;
end

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
% Get middle frame of middle frame cluster
mid_frame_index = fiv(:, round(nB/2));
ref_frame_index = mid_frame_index(round(xB/2));

frame = single(imread(fullfile(octa_path, octa_fname), ...
    ref_frame_index));
frame = subtractBackground(frame);
frame = resampleOCU(frame, p, interpIndex);
fft_frame = ocu_fft(frame);

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

%% Display whole frame for vessel and static tissue ROI selection
frame = applyDispComp(frame, Gc);
fft_frame = ocu_fft(frame);

% Vessel
f = figure;
ax = gca;
imagesc(fft_frame)
title('Select Vessel, Double click roi when done');
vessel_roi = imrect(ax, ...
    [size(frame,2)/3, size(frame,1)/2/3, ...
    size(frame,2)/3, size(frame,1)/2/3]);
vessel_roi = round(wait(vessel_roi));

% Tissue
title('Select nearby static tissue, Double click roi when done');
tissue_roi = imrect(ax, ...
    [size(frame,2)/3, size(frame,1)/2/3, ...
    size(frame,2)/3, size(frame,1)/2/3]);
tissue_roi = round(wait(tissue_roi));
close(f);

%% Get background vector from whole volume
bg = getBG(fullfile(octa_path, octa_fname), wb);

%% Read test frames
frames = zeros(ht, wd, size(fiv, 1), 'single');
for tt=1:xB
    frames(:,:,tt) = single(...
        imread(fullfile(octa_path, octa_fname), ...
        mid_frame_index(tt)));
end

% Process
frames = subtractBackground(frames, bg);
frames = resampleOCU(frames, p, interpIndex);
frames = applyDispComp(frames, Gc);

%% Measure decorrelation ratio between vessel and tissue
d_vessel_to_tissue = zeros(numel(NSPLIT), numel(NBW));
time_mat = d_vessel_to_tissue;
df = figure;
dax = gca;
for ii=1:numel(NSPLIT)
    for jj=1:numel(NBW)
        tic
        %% Set up SSADA filters
        GM = split_spectrum(1:2048, NSPLIT(ii), NBW(jj), 0);
        nGM = size(GM, 1);
        
        % Filter the spectrum
        amp_frames  = zeros(ht/2, wd, xB, nGM, 'single');
        for tt=1:xB
            for mm=1:nGM
%                 filt_frame = single(complex(zeros(ht, wd)));
%                 for xx=1:size(frames, 2)
%                     filt_frame(:,xx) = conv(frames(:,xx,tt), GM(mm, :)', 'same');
%                 end
                filt_frame = frames(:,:,tt) .* GM(mm ,:)';
                fft_frame = fft(filt_frame, [] , 1);
                fft_frame = fft_frame(1:ht/2, :);
                
%                 ang_frame = angle(fft_frame);
%                 unwrapped_frame = unwrap(ang_frame,[],1);
%                 amplitude = abs(fft_frame./exp(1i.*unwrapped_frame));
                amp_frames(:,:,tt,mm) = abs(fft_frame)./2048;
            end
        end

        %% Register the unfiltered amp frames
        % Get unfiltered amp frames
        crop_roi = {25:ht/4, 15:wd}; % todo: figure out a better way
        reg_frames = zeros(...
            numel(crop_roi{1}), ...
            numel(crop_roi{2}), xB, 'single');
        for tt=1:xB
            fft_frame = abs(fft(frames(:,:,tt), [], 1))./ht;
            reg_frames(:,:,tt) = fft_frame(crop_roi{1}, crop_roi{2});
        end
        % Choose reference frame
        % for now, this is just the middle frame
        rfi = round(xB/2); % reference frame index
        ref_frame = reg_frames(:,:,rfi);

        % Register amp frames
        reg_amp_frames = amp_frames;
        fprintf('Frame cluster %i\n', jj);
        for tt=1:xB
            if tt==rfi
                % Don't register the reference frame to itself
                continue;
            end
            mov_frame = reg_frames(:,:,tt);
            movingRefObj = imref2d(size(mov_frame));
            fixedRefObj = imref2d(size(ref_frame));

            % Phase correlation
            tform = imregcorr(mov_frame, movingRefObj, ...
                ref_frame, fixedRefObj, ...
                'transformtype', 'translation', 'Window', true);
            fprintf('dx: %0.2f, dy: %0.2f\n', tform.T(3,1), tform.T(3,1));

            % Reset spatial reference objects to filtered amp frames
            fixedRefObj  = imref2d(size(amp_frames(:,:,1,1)));
            movingRefObj = fixedRefObj;
            for mm=1:nGM
                % Apply transformation to filtered amp frames
                reg_amp_frames(:,:,tt,mm) = imwarp(...
                    amp_frames(:,:,tt,mm), movingRefObj, ...
                    tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
            end
        end
    
        % Measure decorrelation on registered and filtered amp frames
        d = zeros(ht/2, wd, 'single');
        for tt=1:xB-1
            for mm=1:nGM
                d = d + decorrelation(...
                    reg_amp_frames(:,:,tt,mm), ...
                    reg_amp_frames(:,:,tt+1,mm));
            end
        end
        octa_frame = 1 - d./nGM./(xB-1);
        
        % Measure decorrelation in ROIs
        d_vessel = mean(octa_frame(...
            vessel_roi(2):vessel_roi(2)+vessel_roi(4)-1, ...
            vessel_roi(1):vessel_roi(1)+vessel_roi(3)-1), 'all');
        d_tissue = mean(octa_frame(...
            tissue_roi(2):tissue_roi(2)+tissue_roi(4)-1, ...
            tissue_roi(1):tissue_roi(1)+tissue_roi(3)-1), 'all');
        
        d_vessel_to_tissue(ii,jj) = ...
            (d_vessel - d_tissue) / (d_vessel + d_tissue);
        
        time_mat(ii,jj) = toc;
        
        figure(df);
        imagesc(octa_frame);
        imrect(dax, vessel_roi);
        imrect(dax, tissue_roi);
        title(sprintf('D coeff: %1.3f', d_vessel_to_tissue(ii,jj)));
        drawnow();
%         % VISTA x2 (double TR)
%         xTR = 2; % repetition time multiplier
%         nComps = ceil(xB/xTR); % # comparisons
%         vix = zeros(nComps, 2); % VISTA index
%         for xx=1:nComps
%             vix(xx, :) = [xx, xx+xTR];
%         end
% 
%         d = zeros(ht/2, wd, 'single');
%         for tt=1:nComps
%             for mm=1:nGM
%                 d = d + decorrelation(...
%                     reg_amp_frames(:,:,vix(tt,1),mm), ...
%                     reg_amp_frames(:,:,vix(tt,2),mm));
%             end
%         end
%         octa_frame_tr2 = 1 - d./nGM./(nComps);
% 
%         % Threshold based on amplitude
%         % todo: do this based on the whole volume
%         avg_amp_frame = mean(mean(amp_frames, 4), 3);
%         amp_mask = avg_amp_frame >= ...
%             mean(avg_amp_frame(:)) + 2*std(avg_amp_frame(:));
%         octa_frame = octa_frame.*amp_mask;

        waitbar(ii/numel(NSPLIT), wb, ...
            sprintf('# splits: %i, bw: %1.1f', NSPLIT(ii), NBW(jj)));
    end
end






