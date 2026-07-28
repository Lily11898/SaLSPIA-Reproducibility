function setup_paths()
%SETUP_PATHS  Add the SaLSPIA source folders to the MATLAB path.
%   Run this once per session from anywhere:
%       >> run('/path/to/SaLSPIA_reproducibility_archive/setup_paths.m')
%   or, if the archive root is the current folder:
%       >> setup_paths
%
%   Adds algorithms/, utils/ and experiments/ to the path. Data files in
%   data/ are resolved automatically by salspia_data.m, so data/ does not
%   need to be on the path.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'algorithms'));
addpath(fullfile(here, 'utils'));
addpath(fullfile(here, 'experiments'));

fprintf('SaLSPIA paths added:\n  %s\n  %s\n  %s\n', ...
    fullfile(here, 'algorithms'), ...
    fullfile(here, 'utils'), ...
    fullfile(here, 'experiments'));
end
