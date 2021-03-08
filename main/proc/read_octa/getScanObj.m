function [scanObj, fail, err] = getScanObj(octa_ffname)
%getXB finds the number of repeat B-scans used for this scan

%% Defaults
scanObj = [];
fail    = false;
err     = [];

%% Constants
param_names = {...
    'A-scans',              'A';
    'B-scans',              'B';
    'C-scans',              'C';
    'Repeat B-scans',       'xB';
    'Scan Angle (degrees)', 'theta';
    'X FOV (degrees)',      'X';
    'Y FOV (degrees)',      'Y';
    'X offset (degrees)',   'dX';
    'Y offset (degrees)',   'dY';
    'Bidirectional',        'bd'};

    
%% Find scan parameter file
[octa_path, octa_name, ~] = fileparts(octa_ffname);
scan_fname = [octa_name, '-scanobj.txt'];
if exist(fullfile(octa_path, scan_fname), 'file') == 0
    warning('Scan parameters not found, skipping');
    fail = true;
    err = 'file not found';
    return;
end

%% Read file
fid = fopen(fullfile(octa_path, scan_fname), 'r');
try
f_line = '';
while ~isnumeric(f_line)
    % Get parameter name
    if contains(f_line, ':')
        line_parts = strtrim(strsplit(f_line, ':'));
        param_index = strcmpi(param_names(:,1), line_parts{1});
        % Handle booleans
        if isnan(str2double(line_parts{2}))
            line_parts{2} = num2str(eval(lower(line_parts{2})));
        end
        % Convert strings to #'s and add as fields of the scan param object
        eval(sprintf('scanObj.%s = %d;', ...
            param_names{param_index, 2}, ...
            str2double(line_parts{2})));
    end
    f_line = fgetl(fid);
end
catch MException
    fclose(fid);
    err = MException.message;
    fail = true;
    return;
end
fclose(fid);



end

