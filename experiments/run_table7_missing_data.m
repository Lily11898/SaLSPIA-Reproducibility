%% ========================================================================
%  run_table7_missing_data.m
%  Unified reproduction entry for manuscript Table 7.
%
%  Runs:
%    1. the three rank-deficient missing-curve cases (Ex. 4.1, 4.2, 4.7);
%    2. the rank-deficient Franke missing-surface case (Ex. 4.8);
%    3. combines their five-run CPU mean/sample-SD and held-out errors into
%       one CSV and one LaTeX table.
%
%  Outputs:
%    results/table7_missing_data/table7_missing_data.csv
%    results/table7_missing_data/table7_missing_data.tex
%  ========================================================================
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
repo_dir   = fileparts(script_dir);
addpath(fullfile(repo_dir, 'algorithms'));
addpath(fullfile(repo_dir, 'utils'));

fprintf('\n============================================================\n');
fprintf('Unified reproduction of manuscript Table 7\n');
fprintf('============================================================\n');

% This driver produces the three held-out missing-curve cases with five-run
% CPU mean/sample SD. It also regenerates Table 6, which is inexpensive and
% shares the same representative-residual calculations.
run(fullfile(script_dir, 'run_table6_geometric_residuals.m'));

% Recompute paths because the called script begins with "clear".
script_dir = fileparts(mfilename('fullpath'));

% Formal 101x101 Franke experiment with 41x41 controls, 8627 observations,
% and rank(A_obs)=1639 under MATLAB's default SVD tolerance.
run(fullfile(script_dir, 'run_table7_missing_surface.m'));

% Recompute paths again after the called script's "clear".
script_dir = fileparts(mfilename('fullpath'));
repo_dir   = fileparts(script_dir);
out_dir    = fullfile(repo_dir, 'results', 'table7_missing_data');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

curve_file = fullfile(repo_dir, 'results', 'geometric_error_tables', ...
    'missing_heldout_errors.csv');
surface_file = fullfile(repo_dir, 'results', 'missing_data_surface', ...
    'missing_surface_summary.csv');

Tc = readtable(curve_file);
Ts = readtable(surface_file);

vars = {'Example','Type','Method','ObservedPoints','RemovedPoints', ...
    'ControlPoints','Rank','Converged','IT','CPU_mean','CPU_std', ...
    'E_inf','HoleRMS','HoleMax'};
rows = cell(0, numel(vars));

curve_meta = {
    'Ex.~4.1', 267, 64, 101, 96;
    'Ex.~4.2', 458, 63, 201, 189;
    'Ex.~4.7', 206, 63, 101, 90
};
method_order = {'LSPIA','LSPIA-Lin2018','ALSPIA','SaLSPIA'};

for ie = 1:size(curve_meta, 1)
    ex = curve_meta{ie, 1};
    for im = 1:numel(method_order)
        method = method_order{im};
        if strcmp(method, 'LSPIA')
            rows(end+1, :) = {ex, 'curve', method, curve_meta{ie,2}, ...
                curve_meta{ie,3}, curve_meta{ie,4}, curve_meta{ie,5}, ...
                false, 10000, NaN, NaN, NaN, NaN, NaN}; %#ok<SAGROW>
        else
            idx = strcmp(Tc.Example, ex) & strcmp(Tc.Method, method);
            assert(nnz(idx) == 1, 'Expected one curve row for %s / %s.', ex, method);
            rows(end+1, :) = {ex, 'curve', method, Tc.ObservedPoints(idx), ...
                Tc.RemovedPoints(idx), Tc.ControlPoints(idx), Tc.Rank(idx), ...
                true, Tc.IT(idx), Tc.CPU_mean(idx), Tc.CPU_std(idx), ...
                Tc.E_inf(idx), Tc.HoleRMS(idx), Tc.HoleMax(idx)}; %#ok<SAGROW>
        end
    end
end

for im = 1:numel(method_order)
    method = method_order{im};
    idx = strcmp(Ts.name, method);
    assert(nnz(idx) == 1, 'Expected one Franke row for %s.', method);
    converged = Ts.IT(idx) < 10000;
    if converged
        cpu_mean = Ts.CPU_mean(idx);
        cpu_std = Ts.CPU_std(idx);
        e_inf = Ts.Einf(idx);
        hole_rms = Ts.HoleRMS(idx);
        hole_max = Ts.HoleMax(idx);
    else
        cpu_mean = NaN;
        cpu_std = NaN;
        e_inf = NaN;
        hole_rms = NaN;
        hole_max = NaN;
    end
    rows(end+1, :) = {'Ex.~4.8', 'surface', method, 8627, 1574, ...
        1681, 1639, converged, Ts.IT(idx), cpu_mean, cpu_std, ...
        e_inf, hole_rms, hole_max}; %#ok<SAGROW>
end

T = cell2table(rows, 'VariableNames', vars);
csv_file = fullfile(out_dir, 'table7_missing_data.csv');
tex_file = fullfile(out_dir, 'table7_missing_data.tex');
writetable(T, csv_file);
write_table7_latex(T, tex_file);

fprintf('\nComplete Table 7 outputs saved to:\n');
fprintf('  %s\n', csv_file);
fprintf('  %s\n', tex_file);

function write_table7_latex(T, fname)
fid = fopen(fname, 'w');
if fid == -1
    error('Cannot write LaTeX table: %s', fname);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '\\begin{table}[H]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\scriptsize\n');
fprintf(fid, ['\\caption{Rank-deficient missing-data fitting. CPU time (s) is ' ...
    'reported as the mean $\\pm$ sample standard deviation over five runs.}\n']);
fprintf(fid, '\\label{tab:rank_deficient}\n');
fprintf(fid, '\\begin{tabular}{llccccc}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Example & Method & IT & CPU (s), mean $\\pm$ SD & $E_\\infty$ & Hole RMS & Hole Max \\\\\n');
fprintf(fid, '\\hline\n');

for i = 1:height(T)
    if T.Converged(i)
        it_text = sprintf('%d', T.IT(i));
        cpu_text = sprintf('$%.6f \\pm %.6f$', T.CPU_mean(i), T.CPU_std(i));
        e_text = sprintf('$%.2e$', T.E_inf(i));
        rms_text = sprintf('$%.3e$', T.HoleRMS(i));
        max_text = sprintf('$%.3e$', T.HoleMax(i));
    else
        it_text = '$>10^4$';
        cpu_text = '\#';
        e_text = '\#';
        rms_text = '--';
        max_text = '--';
    end
    fprintf(fid, '%s & %s & %s & %s & %s & %s & %s \\\\\n', ...
        T.Example{i}, T.Method{i}, it_text, cpu_text, e_text, ...
        rms_text, max_text);
    if i < height(T) && ~strcmp(T.Example{i}, T.Example{i+1})
        fprintf(fid, '\\hline\n');
    end
end

fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\par\\smallskip\n');
fprintf(fid, ['\\parbox{0.95\\textwidth}{\\footnotesize \\#: failed to converge ' ...
    'within $10^4$ iterations; --: held-out error is not reported for the ' ...
    'non-convergent method.}\n']);
fprintf(fid, '\\end{table}\n');
end
