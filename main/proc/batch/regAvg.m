%% Imports
addpath(genpath('F:\OCT-software\animal-octa\main\proc'));

% Reg avg
[tiff_fnames, tiff_path] = uigetfile('*.tiff', 'Select image', 'F:\img', ...
    'multiselect', 'on');
if ~iscell(tiff_fnames)
    tiff_fnames = {tiff_fnames};
end
tiff_fnames = tiff_fnames';

% All-important waitbar
wbM = waitbar(0, sprintf('Processing %s%s...', tiff_fnames{1}));
wbM.Children.Title.Interpreter = 'none';
waitbar(0, wbM, sprintf('Processing %s...', tiff_fnames{1}));
% Make a within loop waitbar as well
wbS = waitbar(0);

for ii=1:numel(tiff_fnames)
    waitbar(ii/numel(tiff_fnames), wbM, ...
        sprintf('Processing %s...', tiff_fnames{ii}));
    
    tiff_fname = tiff_fnames{ii};
    octa_fname = strrep(tiff_fname, '-fft.tiff', '.octa');
    [xB, fail] = getXB(fullfile(tiff_path, octa_fname));
    if fail || xB < 2
        continue;
    end
    
    % Prep output
    out_fname = strrep(tiff_fname, '.tiff', '-ra.tiff');
    
    % Get indexing info
    tiff_info = imfinfo(fullfile(tiff_path, tiff_fname));
    fiv = 1:numel(tiff_info); % Frame indexing vector
    fiv = reshape(fiv, xB, numel(tiff_info)/xB);
    nB = size(fiv, 2);
    ht = tiff_info(1).Height;
    wd = tiff_info(2).Width;
    
    %% Read the scans as clusters
    for jj=1:nB
        frames = zeros(ht, wd, size(fiv, 1), 'single');
        for tt=1:xB
            frames(:,:,tt) = single(...
                imread(fullfile(tiff_path, tiff_fname), ...
                fiv(tt, jj)));    
        end
        frames = (10.^(frames./65535))-1;
        
        %% Register this cluster of frames
        crop_roi = {25:ht/3, 15:wd}; % todo: figure out a better way
        crop_frames = zeros(...
            numel(crop_roi{1}), ...
            numel(crop_roi{2}), xB, 'single');
        for tt=1:xB
            crop_frames(:,:,tt) = frames(crop_roi{1}, crop_roi{2}, tt);
        end
        % Choose reference frame
        % for now, this is just the middle frame
        rfi = round(xB/2); % reference frame index
        ref_frame = crop_frames(:,:,rfi);
        
        % Register cropped frames
        reg_frames = frames;
        fprintf('Frame cluster %i\n', jj);
        for tt=1:xB
            if tt==rfi
                % Don't register the reference frame to itself
                continue;
            end
            mov_frame = crop_frames(:,:,tt);
            movingRefObj = imref2d(size(mov_frame));
            fixedRefObj = imref2d(size(ref_frame));
            
            % Phase correlation
            tform = imregcorr(mov_frame, movingRefObj, ...
                ref_frame, fixedRefObj, ...
                'transformtype', 'translation', 'Window', true);
            fprintf('dx: %0.2f, dy: %0.2f\n', tform.T(3,1), tform.T(3,1));
            
            % Reset spatial reference objects to original frames
            fixedRefObj  = imref2d(size(frames(:,:,1)));
            movingRefObj = fixedRefObj;
            
            % Apply transformation to original frames
            reg_frames(:,:,tt) = imwarp(...
                frames(:,:,tt), movingRefObj, ...
                tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
        end
        mean_frame = mean(reg_frames, 3);
        mean_frame = uint16(log10(mean_frame+1).*65535);
        
        if jj==1
            wm = 'overwrite';
        else
            wm = 'append';
        end
        imwrite(mean_frame, fullfile(tiff_path, out_fname), ...
            'writemode', wm);
        
        waitbar(jj/nB, wbS);
    end
end

