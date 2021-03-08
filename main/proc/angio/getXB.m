function [xB, fail] = getXB(octa_ffname)
%getXB finds the number of repeat B-scans used for this scan

%% Defaults
xB = 0;
fail = false;

%% Constants
param_name = 'Repeat B-scans:';

%% Find scan parameter file
[octa_path, octa_name, ~] = fileparts(octa_ffname);
scan_fname = [octa_name, '-scanobj.txt'];
if exist(fullfile(octa_path, scan_fname), 'file') == 0
    warning('Scan parameters not found, skipping');
    fail = true;
    return;
end

%% Read file
fid = fopen(fullfile(octa_path, scan_fname), 'r');
f_line = '';
while ~isnumeric(f_line)
    if contains(f_line, param_name, 'ignorecase', true)
        line_parts = strtrim(strsplit(f_line, ':'));
        xB = str2double(line_parts{end});
    end
    f_line = fgetl(fid);
end
fclose(fid);

if xB == 0
    warning('Parameter not found');
    fail = true;
    return;
end


end

