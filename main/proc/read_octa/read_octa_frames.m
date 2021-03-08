function frames = read_octa_frames(octa_ffname, scan, frame_indices, on_gpu)
%todo: document

%% Optional input
if exist('on_gpu', 'var') == 0 || isempty(on_gpu)
    on_gpu = false;
end

%% Get dimensions from scan object
dims = [2048, scan.A*(scan.bd+1), scan.B*scan.xB*scan.C];

%% If frame indices are not given, read whole volume
if exist('frame_indices', 'var') == 0 || isempty(frame_indices)
    frame_indices = 1:dims(3);
end

% Check for oob before op
if any(frame_indices > dims(3))
    error('Frame index out of bounds, check that scan header matches');
end

%% Read frames
%% Read file
% todo: make precision and machine format optional input args
fid  = fopen(octa_ffname, 'r');
frames = zeros(dims(1), dims(2), numel(frame_indices), 'uint16');

%% Transfer to GPU
if on_gpu
    frames = gpuArray(frames);
end

for ii=1:numel(frame_indices)
    % Move to start position of this frame
    status = fseek(fid, ...
        (prod(dims(1:2))*(frame_indices(ii)-1))*2, ... 
        'bof'); % x2 16-bit = 2bytes
    if status ~= 0
        error('I don''t know binary as well as I thought');
    end
    
    % Read frame
    frames(:,:,ii) = fread(fid, [dims(1), dims(2)], 'uint16', 0, 'l');
end

% Close file and have a good day
fclose(fid);

end

