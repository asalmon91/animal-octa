%% Import
addpath(genpath('F:\OCT-software\animal-octa\main\proc'));

%% Get dataset
octa_path = 'F:\img\2019.11.03-DM_180402\OCTA\2019_11_03_OS\Raw';
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

%% Get background vector
bg = single(mean(imread(fullfile(bg_path, bg_fname)), 2));

for ii=1:numel(octa_dir)
    tiff_out_name = strrep(octa_dir(ii).name, '.tiff', '-fft.tiff');
    if exist(fullfile(octa_path, tiff_out_name), 'file') ~= 0
        continue;
    end
    
    %% Get tiff stack info
    octa_info = imfinfo(fullfile(octa_path, octa_dir(ii).name));
    
    %% Optimize dispersion based on the middle frame
    % TODO: include frame selection then ROI selection
    % TODO: ask whether they want to calculate dispersion for just the
    % first scan and apply those parameters to the rest
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
    
    %% Get background vector from whole volume
    bg = getBG(fullfile(octa_path, octa_dir(ii).name), wb);
    
    %% Process the rest of the scans
    for jj=1:numel(octa_info)
        frame = single(imread(fullfile(octa_path, octa_dir(ii).name), jj));
        
        frame = subtractBackground(frame, bg);
        frame = resampleOCU(frame, p, interpIndex, wb);
        frame = applyDispComp(frame, Gc, wb);
        frame = ocu_fft(frame, wb);
        % Intensity transformation for writing to tiff
        frame = uint16(log10(frame./2048 +1) *65535);
        
        if jj==1
            wm = 'overwrite';
        else
            wm = 'append';
        end
        imwrite(frame, fullfile(octa_path, tiff_out_name), ...
            'writemode', wm);
        
        waitbar(jj/numel(octa_info), wb);
    end
    
    waitbar(ii/numel(octa_dir), wb, ...
        sprintf('Done converting %s.', octa_dir(ii).name));
end
close(wb);






