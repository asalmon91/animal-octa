
img_path = 'E:\test\';
img_fname = 'noise-20190923_215913-OS.octa';
fid = fopen(fullfile(img_path, img_fname, 'r'));
dims = fread(fid, 6, 'uint16', 'l'); % First 6 elements are dimensions
I = uint16(fread(fid, prod(dims(dims~=0)), 'uint16', 0, 'l'));
fclose(fid);
I = reshape(I, flip(dims(dims~=0))');

out_tiff_fname = strrep(img_fname, '.octa', '.tiff');
for ii=1:size(I,3)
    if ii==1
        append_mode = 'overwrite';
    else
        append_mode = 'append';
    end
    
    imwrite(I(:,:,ii), fullfile(img_path, out_tiff_fname), ...
        'writemode', append_mode);
end










