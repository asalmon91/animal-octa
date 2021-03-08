function frames = read_octa_frames(octa_ffname, scan, frame_indices)
%todo

%% Constants
% NDIM = 4; % number of elements encoding frame size

%% Read file
% todo: make precision and machine format optional input args
fid  = fopen(octa_ffname, 'r');

% Get matrix size
% dims = fread(fid, NDIM, 'uint16', 'l'); % First 6 elements are dimensions
% dims = flip(dims(1:2:end));


%% Read frames
dims = [2048, scan.A*(scan.bd+1), scan.B*scan.xB*scan.C];
if any(frame_indices > dims(3))
    error('Frame index out of bounds');
end

frames = zeros(dims(1), dims(2), numel(frame_indices), 'uint16');
for ii=1:numel(frame_indices)
    % Move to start position of this frame
    status = fseek(fid, ...
        (prod(dims(1:2))*(frame_indices(ii)-1))*2, ... 
        'bof'); % x2 16-bit = 2bytes
    if status ~= 0
        error('I don''t know binary as well as I thought');
    end
    
    % Read frame
    frames(:,:,ii) = uint16(...
        fread(fid, [dims(1), dims(2)], 'uint16', 0, 'l'));

end

% Close file and have a good day
fclose(fid);

end

