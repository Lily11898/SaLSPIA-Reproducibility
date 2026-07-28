%% ========================================================================
%  run_table2_parameter_sensitivity_curves.m
%  Reproduce the SaLSPIA curve parameter-sensitivity study in Table 2.
%
%  Cases:
%    1. Curve Ex. 4.3: full-rank reindeer contour fitting
%    2. Missing G-loop: rank-deficient fitting with four missing regions
%
%  Only SaLSPIA is tested. This script reproduces the sigma and delta
%  statements in the discussion and the first two rows of Table 2 for M.
%  ========================================================================
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(script_dir);
addpath(fullfile(repo_dir, 'algorithms'));
addpath(fullfile(repo_dir, 'utils'));

out_dir = fullfile(repo_dir, 'results', 'parameter_sensitivity');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

p_deg   = 3;
tol_Ek  = 1e-6;
maxiter = 10000;

sigma_values = [0.5, 1e-1, 1e-2, 1e-4, 1e-6];
delta_values = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2];
M_values     = [1, 3, 5, 10, 20];

defaults.c           = 1e-4;
defaults.M           = 10;
defaults.delta       = 1e-8;
defaults.eps_saf     = 1e-30;
defaults.p           = 3;
defaults.e_tol       = 1e-15;
defaults.max_halving = 100;
defaults.mode        = 'salspia';

cases = setup_cases(repo_dir, p_deg);
summary_rows = {};
infos = cell(numel(cases), 3);

fprintf('\n============================================================\n');
fprintf('SaLSPIA parameter sensitivity\n');
fprintf('============================================================\n');

for ic = 1:numel(cases)
    cs = cases(ic);
    fprintf('\n%s\n', cs.label);
    fprintf('  data=%d, control=%d, rank=%d\n', ...
        size(cs.Q, 1), size(cs.P0, 1), cs.rankA);

    infos{ic, 1} = cell(numel(sigma_values), 1);
    for iv = 1:numel(sigma_values)
        params = defaults;
        params.c = sigma_values(iv);
        [~, info] = ablation_lspia(cs.A, cs.Q, cs.P0, tol_Ek, maxiter, params);
        infos{ic, 1}{iv} = info;
        summary_rows(end+1, :) = make_row(cs, 'sigma', sigma_values(iv), ...
            params, info); %#ok<SAGROW>
        print_row('sigma', sigma_values(iv), info);
    end

    infos{ic, 2} = cell(numel(delta_values), 1);
    for iv = 1:numel(delta_values)
        params = defaults;
        params.delta = delta_values(iv);
        [~, info] = ablation_lspia(cs.A, cs.Q, cs.P0, tol_Ek, maxiter, params);
        infos{ic, 2}{iv} = info;
        summary_rows(end+1, :) = make_row(cs, 'delta', delta_values(iv), ...
            params, info); %#ok<SAGROW>
        print_row('delta', delta_values(iv), info);
    end

    infos{ic, 3} = cell(numel(M_values), 1);
    for iv = 1:numel(M_values)
        params = defaults;
        params.M = M_values(iv);
        [~, info] = ablation_lspia(cs.A, cs.Q, cs.P0, tol_Ek, maxiter, params);
        infos{ic, 3}{iv} = info;
        summary_rows(end+1, :) = make_row(cs, 'M', M_values(iv), ...
            params, info); %#ok<SAGROW>
        print_row('M', M_values(iv), info);
    end
end

T = cell2table(summary_rows, 'VariableNames', ...
    {'Case', 'Sweep', 'Value', 'DataPoints', 'ControlPoints', 'Rank', ...
     'Sigma', 'M', 'Delta', 'IT', 'CPU', 'E_inf', 'FittingError', ...
     'TotalHalving', 'Status'});

writetable(T, fullfile(out_dir, 'parameter_sensitivity_summary.csv'));
save(fullfile(out_dir, 'parameter_sensitivity_results.mat'), ...
    'T', 'infos', 'cases', 'sigma_values', 'delta_values', 'M_values', ...
    'defaults', 'tol_Ek', 'maxiter');
write_table2_latex(T, fullfile(out_dir, 'table2_parameter_sensitivity.tex'));

disp(T);
fprintf('\nSaved outputs to:\n  %s\n', out_dir);

%% ===== Case setup =======================================================
function cases = setup_cases(repo_dir, p_deg)
cases = struct('label', {}, 'short_name', {}, 'Q', {}, 'A', {}, ...
    'P0', {}, 'rankA', {});

% Curve Ex. 4.3: resample the 508-point source contour to 1000 points.
n = 287;
m1 = 1000;
Q = load(salspia_data('cur_data deer'));
if size(Q, 1) ~= m1
    idx = round(linspace(1, size(Q, 1), m1));
    Q = Q(idx, :);
end
[t, knots] = make_params_knots(Q, n, p_deg);
A = build_collocation(knots, p_deg, t);
P0 = select_initial_ctrl_pts(Q, n);
cases(1) = make_case('Curve Ex. 4.3', 'Curve_Ex4_3', Q, A, P0);

% Missing G-loop: same manuscript-faithful construction as Tables 4 and 9.
n = 100;
Q_full = load(salspia_data('s_loop_curve_data.txt'));
[t_full, knots] = make_params_knots(Q_full, n, p_deg);
[Q, t_keep] = remove_by_param_intervals( ...
    Q_full, t_full, [0.15, 0.38, 0.62, 0.99], 0.03);
A = build_collocation(knots, p_deg, t_keep);
P0 = select_initial_ctrl_pts(Q_full, n);
cases(2) = make_case('Missing G-loop', 'Missing_Gloop', Q, A, P0);
end

function cs = make_case(label, short_name, Q, A, P0)
cs.label = label;
cs.short_name = short_name;
cs.Q = Q;
cs.A = A;
cs.P0 = P0;
cs.rankA = rank(A);
end

function [Q_keep, t_keep] = remove_by_param_intervals( ...
    Q_full, t_full, hole_centers, hole_half_width)
keep = true(size(t_full));
for ih = 1:numel(hole_centers)
    lo = hole_centers(ih) - hole_half_width;
    hi = hole_centers(ih) + hole_half_width;
    keep((t_full >= lo) & (t_full <= hi)) = false;
end
Q_keep = Q_full(keep, :);
t_keep = t_full(keep);
end

%% ===== Output helpers ===================================================
function row = make_row(cs, sweep, value, params, info)
row = {cs.short_name, sweep, value, size(cs.Q, 1), size(cs.P0, 1), ...
    cs.rankA, params.c, params.M, params.delta, info.iter_count, ...
    info.cpu_time, info.Ek_history(end), info.err_history(end), ...
    info.total_halving, info.status};
end

function print_row(sweep, value, info)
fprintf('  %-5s=%-8.1e IT=%-3d E=%.2e halv=%-3d status=%s\n', ...
    sweep, value, info.iter_count, info.Ek_history(end), ...
    info.total_halving, info.status);
end

function write_table2_latex(T, fname)
curve = T(strcmp(T.Case, 'Curve_Ex4_3') & strcmp(T.Sweep, 'M'), :);
gloop = T(strcmp(T.Case, 'Missing_Gloop') & strcmp(T.Sweep, 'M'), :);

fid = fopen(fname, 'w');
if fid == -1
    warning('Cannot write LaTeX table: %s', fname);
    return;
end
fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Sensitivity to the memory length $M$ with $\\sigma=10^{-4}$ and $\\delta=10^{-8}$ fixed.}\n');
fprintf(fid, '\\begin{tabular}{llrrrrr}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Example & Metric & $M=1$ & $M=3$ & $M=5$ & $M=10$ & $M=20$ \\\\\n');
fprintf(fid, '\\hline\n');
write_case_rows(fid, 'Curve Ex. 4.3', curve);
write_case_rows(fid, 'Missing G-loop', gloop);
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
end

function write_case_rows(fid, name, T)
fprintf(fid, '%s & $E_\\infty$ ', name);
fprintf(fid, '& %.2e ', T.E_inf);
fprintf(fid, '\\\\\n');
fprintf(fid, ' & IT ');
fprintf(fid, '& %d ', T.IT);
fprintf(fid, '\\\\\n');
fprintf(fid, ' & Halvings ');
fprintf(fid, '& %d ', T.TotalHalving);
fprintf(fid, '\\\\\n');
end
