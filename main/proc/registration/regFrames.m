function reg_frames = regFrames(frames, reg_method, tol)
%regFrames Registers OCT frames
%   First uses a custom optimization approach; optimizes translation and
%   vertical shear to minimize image MSE. If that fails,
%   phase-correlation-based translation is used.

%% Optional input
if exist('reg_method', 'var') == 0 || isempty(reg_method)
    reg_method = 'fast';
else
    if ~strcmpi(reg_method, 'fast') && ~strcmpi(reg_method, 'accurate')
        error('Input: %s not recognized, options are "fast" and "accurate"', reg_method);
    end
end

%% Check required input
% fminsearch is not compatible with gpu-arrays
if ~strcmpi(reg_method, 'fast') && isa(frames, 'gpuArray')
    frames = gather(frames);
end
[~, ~, xRpt] = size(frames);

%% Optional inputs
if strcmpi(reg_method, 'accurate') && ...
        (exist('tol', 'var') == 0 || isempty(tol))
    tol = 1e-6;
end

%% Register the linear amplitude frames
reg_frames = frames;

%% Choose reference frame
% for now, this is just the middle frame
rfi = round(xRpt/2); % reference frame index
% ref_frame = reg_roi_frames(:,:,rfi);
ref_frame = frames(:,:,rfi);
dx_dy_shy = [0,0,0];
% Register amp frames
for tt=1:xRpt
    if tt==rfi
        % Don't register the reference frame to itself
        continue;
    end
    mov_frame = frames(:,:,tt);
    
    %% Shear-based registration optimization
    if strcmpi(reg_method, 'accurate')
        [reg_frames(:,:,tt), dx_dy_shy, success] = ...
            reg_OCT_fmin(ref_frame, mov_frame, dx_dy_shy, tol);
        if ~success
            warning('OCT registration failed, trying an alternate method');
        end
    end
    
    %% Phase-correlation-based translation registration
    if strcmpi(reg_method, 'fast') || ~success
        fixedRefObj     = imref2d(size(ref_frame));
        movingRefObj    = imref2d(size(mov_frame));
        % Phase correlation
        warn_msg = ''; %#ok<NASGU>
        tform = imregcorr(mov_frame, movingRefObj, ...
            ref_frame, fixedRefObj, ...
            'transformtype', 'translation', 'Window', true);
        warn_msg = lastwarn();
        if ~strcmp(warn_msg, '')
            tform.T(3, 1:2) = 0;
            %warning('Phase correlation failed, not registering');
        else
            reg_frames(:,:,tt) = imwarp(...
                mov_frame, movingRefObj, ...
                tform, 'OutputView', fixedRefObj, ...
                'SmoothEdges', true, 'interp', 'cubic', ...
                'fillvalues', mean(mov_frame(:)));
        end
    end
end
reg_frames(reg_frames < 0) = 0; % Can exceed limit due to cubic interp

end




