function C_vec = dispComp_fminsearch(img, C2, C3)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Spectrometer calibration
new_index = loadSpecCal();

%% Starting points for dispersion compensation coefficients
% (if not given)
if exist('C2', 'var') == 0
    C2 = -3e-5;
end
if exist('C3', 'var') == 0
    C3 = 1e-9;
end

%% Preallocate FFT matrix
IdFa = zeros(size(img,1)/2, size(img,2), 'double');

%% Get background for subtraction
xm = (mean(double(img),2));

%% Preallocate vector for sharpness values
sharps = zeros(size(img, 2), 1, 'double');

%% fminbnd ranges
C_lims = [-5e-5, -1e-5; 9e-10, 2e-9]; % Empirically determined

%% Get sharpness with default values
sharp_0 = opt_DispComp(C2,img,new_index, IdFa,xm,sharps,true,C3);
sharp_n1 = opt_DispComp(C_lims(1,1), ...
    img, new_index, IdFa,xm,sharps,true, C_lims(2,1));
sharp_p1 = opt_DispComp(C_lims(1,2), ...
    img, new_index, IdFa,xm,sharps,true, C_lims(2,2));
% Not sure if this is how it should be done... The tolerance is about how
% much it varies... Maybe it would be better to try 3 values of C2 and
% base the tolerance on that
tol_fun = 10^(floor(log10(range([sharp_n1, sharp_0, sharp_p1]))));

%% Start with the default optimization options
options = optimset;

%% Modify options setting for C2
% Testing fminsearch vs fminbnd
search_times = [0,0];
bound_times = [0,0];


%% Set up optimization
C_vec = [C2, C3];
C_search_vec = [0,0];
C_bound_vec = [0,0];
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
    
    tic
    % Run optimization
    [C_vec(ii), fval, exitflag, output] = fminsearch(fun, CX0, options);
    search_times(ii) = toc;
    C_search_vec(ii) = C_vec(ii);
    
    tic
    [C_out, fval, exitflag, output] = fminbnd(fun, ...
        C_lims(ii,1), C_lims(ii,2), options);
    bound_times(ii) = toc;
    C_bound_vec(ii) = C_out;
end

figure;
subplot(4,2,1);
bar(categorical({'fminsearch', 'fminbnd'}), ...
    [C_search_vec(1), C_bound_vec(1)]);
ylabel('C2 estimate');

subplot(4,2,2);
bar(categorical({'fminsearch', 'fminbnd'}), ...
    [C_search_vec(2), C_bound_vec(2)])
ylabel('C3 estimate');

subplot(4,2,3);
bar(categorical({'fminsearch', 'fminbnd'}), ...
    [search_times(1), bound_times(1)])
ylabel('C2 time (s)');

subplot(4,2,4);
bar(categorical({'fminsearch', 'fminbnd'}), ...
    [search_times(2), bound_times(2)])
ylabel('C3 time (s)');

%% Generate images based on these results
[~, search_im_out] = opt_DispComp(C_search_vec(1), ...
    img, new_index, IdFa, xm, sharps, true, C_search_vec(2));

[~, bound_im_out] = opt_DispComp(C_bound_vec(1), ...
    img, new_index, IdFa, xm, sharps, true, C_bound_vec(2));

% Manual roi
roi_x = 169;

subplot(4,2,[5,7]);
imagesc(search_im_out);
title('fminsearch');
hold on;
plot([roi_x, roi_x], [1, 1024], '-r');
hold off;

subplot(4,2,[6,8]);
imagesc(bound_im_out);
title('fminbnd');
hold on;
plot([roi_x, roi_x], [1, 1024], '-r');
hold off;

figure;
hold on;
plot(bound_im_out(:, roi_x));
plot(search_im_out(:, roi_x));
hold off;
legend({'fminbnd', 'fminsearch'});
axis tight
set(gca, 'tickdir', 'out');
xlabel('Row');
ylabel('Intensity');

a = [sum(bound_im_out(:, roi_x)), sum(search_im_out(:, roi_x))]

