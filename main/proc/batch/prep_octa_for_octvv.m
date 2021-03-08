%% Constants
crop = [30, 24, 2];

%% Prep OCTA output for OCTVV
in_path = uigetdir('.', 'Select Processed directory');
ra_dir = dir(fullfile(in_path, '*-ra.tiff'));
var_dir = dir(fullfile(in_path, '*-var.tiff'));
ada_dir = dir(fullfile(in_path, '*-ada1.tiff'));

%% Process registered & averaged volumes
parfor ii=1:numel(ra_dir)
    tiff_fname = ra_dir(ii).name;
    [~, tiff_name, tiff_ext] = fileparts(fullfile(in_path, tiff_fname));
    
    % Prep output
    out_fname = [tiff_name, '-jpg.avi'];
    if exist(fullfile(in_path, out_fname), 'file') ~= 0
        warning('%s exists, skipping', out_fname);
        continue;
    end
    
    % Get volume data
    tiff_data = imfinfo(fullfile(in_path, tiff_fname));
    
    % Prep volume
    avi_out = zeros(...
        tiff_data(1).Height, ...
        tiff_data(1).Width, ...
        numel(tiff_data), 'uint16');
    
    for jj=1:numel(tiff_data)
        avi_out(:,:,jj) = imread(fullfile(in_path, tiff_fname), jj);
    end
    
    % Crop
    avi_out = avi_out(crop(1):512, crop(2):end, crop(3):end); %#ok<PFBNS>
    
    % Contrast stretch
    cs_range = stretchlim(avi_out(:));
    for jj=1:size(avi_out, 3)
        avi_out(:,:,jj) = imadjust(avi_out(:,:,jj), cs_range);
    end
    avi_out = uint8(single(avi_out)./65535.*255);
    
    vw = VideoWriter(fullfile(in_path, out_fname), 'motion jpeg avi'); %#ok<TNMLP>
    open(vw);
    for jj=1:size(avi_out, 3)
        writeVideo(vw, avi_out(:,:,jj));
    end
    close(vw);
end

% TODO: FIX DRY VIOLATION
%% Process variance volume
parfor ii=1:numel(var_dir)
    tiff_fname = var_dir(ii).name;
    [~, tiff_name, tiff_ext] = fileparts(fullfile(in_path, tiff_fname));
    
    % Prep output
    out_fname = [tiff_name, '-jpg.avi'];
    if exist(fullfile(in_path, out_fname), 'file') ~= 0
        warning('%s exists, skipping', out_fname);
        continue;
    end
    
    % Get volume data
    tiff_data = imfinfo(fullfile(in_path, tiff_fname));
    
    % Prep volume
    avi_out = zeros(...
        tiff_data(1).Height, ...
        tiff_data(1).Width, ...
        numel(tiff_data), 'uint16');
    
    for jj=1:numel(tiff_data)
        avi_out(:,:,jj) = imread(fullfile(in_path, tiff_fname), jj);
    end
    
    % Crop
    avi_out = avi_out(crop(1):512, crop(2):end, crop(3):end); %#ok<PFBNS>
    
    % Contrast stretch
    cs_range = double([median(avi_out(:)), max(avi_out(:))])./65535;
    for jj=1:size(avi_out, 3)
        avi_out(:,:,jj) = imadjust(avi_out(:,:,jj), cs_range);
    end
    avi_out = uint8(single(avi_out)./65535.*255);
    
    vw = VideoWriter(fullfile(in_path, out_fname), 'motion jpeg avi'); %#ok<TNMLP>
    open(vw);
    for jj=1:size(avi_out, 3)
        writeVideo(vw, avi_out(:,:,jj));
    end
    close(vw);
end

%% Process registered & averaged volumes
parfor ii=1:numel(ada_dir)
    tiff_fname = ada_dir(ii).name;
    [~, tiff_name, tiff_ext] = fileparts(fullfile(in_path, tiff_fname));
    
    % Prep output
    out_fname = [tiff_name, '-jpg.avi'];
    if exist(fullfile(in_path, out_fname), 'file') ~= 0
        warning('%s exists, skipping', out_fname);
        continue;
    end
    
    % Get volume data
    tiff_data = imfinfo(fullfile(in_path, tiff_fname));
    
    % Prep volume
    avi_out = zeros(...
        tiff_data(1).Height, ...
        tiff_data(1).Width, ...
        numel(tiff_data), 'uint16');
    
    for jj=1:numel(tiff_data)
        avi_out(:,:,jj) = imread(fullfile(in_path, tiff_fname), jj);
    end
    
    % Crop
    avi_out = avi_out(crop(1):512, crop(2):end, crop(3):end); %#ok<PFBNS>
    
    % Contrast stretch
%     cs_range = stretchlim(avi_out(:));
%     for jj=1:size(avi_out, 3)
%         avi_out(:,:,jj) = imadjust(avi_out(:,:,jj), cs_range);
%     end
    avi_out = uint8(single(avi_out)./65535.*255);
    
    vw = VideoWriter(fullfile(in_path, out_fname), 'motion jpeg avi'); %#ok<TNMLP>
    open(vw);
    for jj=1:size(avi_out, 3)
        writeVideo(vw, avi_out(:,:,jj));
    end
    close(vw);
end




