% export_data = export(HST)  Export method for HST data.
%
% Usage:
%   export_data = HST.export
%
% This method exports data contained in the HSTData class to a standard
% MATLAB structure array format.
%
% Inputs:
%   HST     HSTData class object/array of objects
%
% Outputs:
%   export_data     Structure array of data using standard data types.
%
% Copyright (C) by Alan Degenhart and Erinn Grigsby
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com


function export_data = export(HST)

% Specify list of valid classes.  If a field is not one of these types,
% then call the export method for that class.
valid_class = { ...
    'cell', ...
    'char', ...
    'double', ...
    'int16', ...
    'int32', ...
    'int8', ...
    'logical', ...
    'single', ...
    'struct',...
    'uint32'};

% Get field names
field_names = fieldnames(HST);
n_fields = length(field_names);

% Initialize structure
export_data = struct();
for i = 1:n_fields
    export_data.(field_names{i}) = [];
end
sz = size(HST);
export_data = repmat(export_data, sz);

% Check the number of elements.  Currently, this function supports up to 2D
% arrays of objects, and will generate an error for larget objects.
assert(numel(sz) <= 2, 'Export only supports up to 2D arrays.')

% Loop over elements
for i = 1:sz(1)
    for j = 1:sz(2)
        % Iterate over fields
        for k = 1:n_fields
            fn = field_names{k};
            if ismember(class(HST(i, j).(fn)), valid_class)
                % If the data is of a standard format, add it to the
                % structure
                export_data(i, j).(fn) = HST(i, j).(fn);
            else
                % If the data is of a different format (most likely a
                % custom class, then call the 'export' method
                try
                    export_data(i, j).(fn) = HST(i, j).(fn).export;
                catch
                    warning('Could not export %s.\n', fn)
                end
            end
        end
    end
end
