function C_vec = dispComp_fminbnd(img, C2, C3, roi)
%% Optimizes dispersion compensation
% img should already be background subtracted and interpolated to linear
% k-space

% %% Spectrometer calibration
% new_index = loadSpecCal();

%% Starting points for dispersion compensation coefficients
% (if not given)
if exist('C2', 'var') == 0
    C2 = -3e-5;
end
if exist('C3', 'var') == 0
    C3 = 1e-9;
end
if exist('roi', 'var') == 0
    roi = [1, 1, size(img,2), floor(size(img,1)/2)];
end

%% Preallocate FFT matrix
IdFa = zeros(size(img,1)/2, size(img,2), 'double');

% %% Get background for subtraction
% xm = (mean(double(img),2));

%% Preallocate vector for sharpness values
sharps = zeros(roi(3), 1, 'double');

%% fminbnd ranges
C_lims = [-5e-5, -1e-5; 9e-10, 2e-9]; % Empirically determined

%% Get sharpness with default values
sharp_n1 = opt_DispComp(C_lims(1,1), ...
    img, new_index, IdFa,xm,sharps,true, C_lims(2,1));
sharp_p1 = opt_DispComp(C_lims(1,2), ...
    img, new_index, IdFa,xm,sharps,true, C_lims(2,2));
tol_fun = 10^(floor(log10(range([sharp_n1, sharp_p1]))));

%% Start with the default optimization options
options = optimset;

%% Set up optimization
C_vec = [C2, C3];
for ii=1:numel(C_vec)
    % Get current starting point
    CX0 = C_vec(ii);
    CY0 = C_vec(C_vec ~= CX0);
    
    % Set switch
    C2_switch = ii==1;
    
    % Determine tolerance
    tol_x = 10^floor(log10(abs(CX0)));
    
    % Modify options setting
    options = optimset(options,'TolFun', tol_fun);
    options = optimset(options,'TolX', tol_x);
    
    % Define anonymous function to declare which variable needs to be
    % optimized
    fun = @(x) opt_DispComp(x, ...
        img, new_index, IdFa, xm, sharps, C2_switch, CY0);
    [C_vec(ii), fval, exitflag, output] = fminbnd(fun, ...
        C_lims(ii,1), C_lims(ii,2), options);
end


