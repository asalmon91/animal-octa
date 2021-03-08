[tiff_fnames, tiff_path] = uigetfile('*.tiff', 'Select tiff stack', ...
    'multiselect', 'on');
if ~iscell(tiff_fnames)
    tiff_fnames = {tiff_fnames};
end
tiff_fnames = tiff_fnames';

wb = waitbar(0, sprintf('Reading %s...', tiff_fnames{1}));
wb.Children.Title.Interpreter = 'none';
    
for ii=1:numel(tiff_fnames)
    tiff_fname = tiff_fnames{ii};
    img_info = imfinfo(fullfile(tiff_path, tiff_fname));
    
    crop_roi = {25:img_info(1).Height, 14:img_info(1).Width, ...
        1:numel(img_info)};
    
    oct = zeros(img_info(1).Height, img_info(2).Width, numel(img_info), ...
        'uint16');
    for jj=1:numel(img_info)
        oct(:,:,jj) = imread(fullfile(tiff_path, tiff_fname), jj);
        
        if mod(jj, 10) == 0 || jj == numel(img_info)
            waitbar(jj/numel(img_info), wb);
        end
    end
    % Crop
    oct = oct(crop_roi{1}, crop_roi{2}, crop_roi{3});
    % Linearize
    oct_lin = (10.^(double(oct)./25e3)-1) .*2048;
    
    % Flip reverse scan and average
    % oct_fwd = oct_lin(:,1:size(oct_lin, 2)/2, :);
    % oct_rev = flip(oct_lin(:,size(oct_lin, 2)/2+1:end, :), 2);
    % oct_lin = oct_fwd;
    % for ii=1:size(oct_fwd, 3)
    %     oct_lin(:,:,ii) = mean(cat(3, oct_fwd(:,:,ii), oct_rev(:,:,ii)), 3);
    % end
    
    % Generate summed volume projection
    mvp = squeeze(sum(oct_lin, 1))';
    % Median filter
    mvp = medfilt2(mvp);
    % Adjust range to 16-bit
    mvp = mvp-min(mvp(:));
    mvp = uint16(mvp./max(mvp(:)).*65535);
    % Contrast stretch
    mvp = imadjust(mvp, stretchlim(mvp));
    
    figure;
    imagesc(mvp);
    axis equal tight off
    
    % Write out
    imwrite(mvp, fullfile(tiff_path, ...
        strrep(tiff_fname, '.tiff', 'svp.tif')));
    
end
close(wb);
