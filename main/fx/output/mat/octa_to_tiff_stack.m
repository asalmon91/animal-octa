function octa_to_tiff_stack(octa_ffname)
%octa_to_tiff_stack converts a binary octa file to a tiff_stack for easy
%viewing in ImageJ or similar

[octa_path, octa_name, ~] = fileparts(octa_ffname);

% Read file
% todo: make precision and machine format optional input args
fid = fopen(octa_ffname, 'r');
dims = fread(fid, 6, 'uint16', 'l'); % First 6 elements are dimensions
I = uint16(fread(fid, prod(dims(dims~=0)), 'uint16', 0, 'l'));
fclose(fid);
I = reshape(I, flip(dims(dims~=0))');
% I_t = permute(I_t, [3,2,1]);

out_tiff_fname = [octa_name, '.tiff'];
for ii=1:size(I,3)
    if ii==1
        append_mode = 'overwrite';
    else
        append_mode = 'append';
    end

    imwrite(I(:,:,ii), fullfile(octa_path, out_tiff_fname), ...
        'writemode', append_mode, 'compression', 'none');
end

end

