function p = salspia_data(name)
%SALSPIA_DATA  Absolute path to a file in the repository's data/ folder.
%   p = salspia_data('cur_data deer') returns the full path so that scripts
%   in experiments/ can load data from ../data/ regardless of the current
%   working directory. Run setup_paths.m once so this helper is on the path.
here = fileparts(mfilename('fullpath'));      % .../utils
p    = fullfile(fileparts(here), 'data', name);
end
