%% Imports
addpath('F:\OCT-software\animal-octa\main\fx\output\mat');

%% Get Directory
octa_path = uigetdir('F:\img', 'Select directory');
if isnumeric(octa_path) && octa_path == 0
    return;
end
% octa_path = 'F:\img\2019.10.22-DM_180402\OCTA\2019_10_22_OS\Calibration';
octa_dir = dir(fullfile(octa_path, '*.octa'));

%% Waitbar
wb = waitbar(0, sprintf('Converting %s%s...', ...
    octa_dir(1).name, octa_dir(1).name));
wb.Children.Title.Interpreter = 'none';
waitbar(0, wb, sprintf('Converting %s...', octa_dir(1).name));

%% Convert binary files
for ii=1:numel(octa_dir)
    tiff_out_name = strrep(octa_dir(ii).name, '.octa', '.tiff');
    if exist(fullfile(octa_path, tiff_out_name), 'file') ~= 0
        continue;
    end
    octa_to_tiff_stack(fullfile(octa_path, octa_dir(ii).name));
    
    waitbar(ii/numel(octa_dir), wb, ...
        sprintf('Done converting %s.', octa_dir(ii).name));
end
close(wb);






