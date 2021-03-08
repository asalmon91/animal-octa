function [oct_frame, octa_frame, Gc] = proc_buffered_oct_gpu(...
    in_frames, A, Gc, p, interpIndex)
%proc_oct_gpu Processes an OCT-A volume on the GPU
% @octa_frames is an M*(A*xB)*1 uint16 array of unprocessed camera images
% where M is the length of the line-scan camera, @A is the # of
% A-scans/B-scan, and xB is the # of repeated B-scans
% @A is the number of A-scans/B-scan, if not input, it's assumed that xB is
% 1, and the width of octa_frames is equal to A.
% @Gc is an Mx1 complex double vector for dispersion compensation. It is
% optional.

% If there are repeated B-scans, it outputs a registered averaged 
% structural frame (@oct_frame) and a registered averaged aniogram 
% (@octa_frame)frame; otherwise it outputs the raw processed frame for 
% @oct_frame, and nothing % for @octa_frame

%% Imports
addpath(genpath('.'));

%% Constants
% Registration tolerance, 1e-4 works okay and takes ~0.5s on my system,
% 1e-6 works quite well but takes 6s
REG_TOL = 1e-4; 

%% Check input, get dimensions
if exist('A','var')==0 || isempty(A)
    xRpt = 1;
    A = size(in_frames, 2);
else
    xRpt = size(in_frames, 2) / A;
    if mod(xRpt, 1) ~= 0
        error('Matrix dimensions don''t match. A=%i', A);
    end
end
ht = size(in_frames, 1);

%% Send input frames to the GPU
try
    in_frames = single(gpuArray(in_frames));
catch merr
    % Catch Out-of-memory error, reset device, try again. If that fails,
    % process on CPU
    if ~strcmp(merr.identifier, 'parallel:gpu:array:OOM')
        rethrow(merr);
    end
    % Reset
    gpuDevice(1);
    try
        in_frames = single(gpuArray(in_frames));
    catch
        if ~strcmp(merr.identifier, 'parallel:gpu:array:OOM')
            rethrow(merr);
        end
    end
end    

%% Load spectrometer calibration
z = ht;
if exist('p','var')==0 || isempty(p)
    p = 1:z;
end
if exist('interpIndex','var')==0 || isempty(interpIndex)
    interpIndex = loadSpecCal();
end

%% Get dispersion compensation
if exist('Gc', 'var') == 0 || isempty(Gc)
    Gc = user_assisted_dispComp([], [], p, interpIndex, in_frames);
end

%% DC Mitigation (background subtraction)
in_frames = in_frames - mean(in_frames, 2);

%% Resampling, Dispersion Compensation, and FFT
in_frames = fft(...
    interp1(p', in_frames, interpIndex', 'linear') ... %interp
    .* Gc', ... % dispComp
    [] , 1); % fft
in_frames = abs(in_frames(1:ht/2, :))./ht; % Get amplitude

%% Check if repeated B-scans
if xRpt > 1
    %% Register frames
    in_frames = reshape(in_frames, [ht/2, A, xRpt]);
    reg_frames = regFrames(in_frames, REG_TOL);
    % At this point reg_frames is no longer on the GPU
    % regFrames is not compatible with gpuArrays
    %% Registered-average
    oct_frame = mean(reg_frames, 3);
    
    %% Get angiogram
    octa_frame = uint16(get_fsada(reg_frames).*(2^16-1));
    
else % Just return oct_frame
    if isa(in_frames, 'gpuArray')
        in_frames = gather(in_frames);
    end
    oct_frame   = in_frames;
    octa_frame  = [];
end


end

