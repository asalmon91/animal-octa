function [reg_mov_frame, dx_dy_shy, success] = ...
    reg_OCT_fmin(ref_frame, mov_frame, dx_dy_shy, opt_tol)
%reg_OCT_fmin Determines the optimal registration of two OCT images
%   Uses only translation and vertical shear

%% Defaults
success = false;

%% Constants
if exist('opt_tol', 'var') == 0 || isempty(opt_tol)
    opt_tol  = 1e-6;
end
MAX_DISP    = 10; % pixels
MAX_SHEAR   = 0.05;

%% Determine error range
mse0 = immse(ref_frame, mov_frame);
mse1 = opt_OCT_reg([-MAX_DISP, -MAX_DISP, MAX_SHEAR], ref_frame, mov_frame);
tol_fun = 10^(floor(log10(abs(mse0-mse1)))-1);

%% Set optimization options
options = optimset('TolFun', tol_fun, 'TolX', opt_tol);

fun = @(x) opt_OCT_reg(x, ref_frame, mov_frame);
% tic
[dx_dy_shy,~,exitflag] = fminsearch(fun, dx_dy_shy, options);
% toc
% todo: handle poor registrations
if any(dx_dy_shy > [MAX_DISP, MAX_DISP, MAX_SHEAR]) || ~exitflag
    warning('Registration failed');
    dx_dy_shy = [0,0,0];
    reg_mov_frame = mov_frame;
    return;
else
    success = true;
end

%% Warp moving frame to match reference frame
tform = affine2d();
tform.T(3,1) = dx_dy_shy(1);
tform.T(3,2) = dx_dy_shy(2);
tform.T(1,2) = dx_dy_shy(3);
fixedRefObj = imref2d(size(ref_frame));
moveRefObj  = imref2d(size(mov_frame));
reg_mov_frame = imwarp(mov_frame, moveRefObj, ...
        tform, 'OutputView', fixedRefObj, ...
        'SmoothEdges', true, 'interp', 'cubic', ...
        'fillvalues', mean(mov_frame(:)));

end

