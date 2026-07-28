%% ========================================================================
%  run_table6_geometric_residuals.m
%  Generate the compact geometric-error tables used to answer reviewers:
%
%  Table A: representative RMS / maximum pointwise residuals for SaLSPIA
%           on one full-rank curve example and one full-rank surface example.
%  Table B: held-out RMS / maximum pointwise errors on the removed data
%           points for the missing-data examples.
%
%  Run setup_paths.m once before this script.
%  ========================================================================
clear; clc; close all;

%% ===== Settings =========================================================
p_deg   = 3;
tol_Ek  = 1e-6;
maxiter = 10000;
n_repeats = 5;                 % used only for the compact held-out comparison

salspia_params.c     = 1e-4;
salspia_params.M     = 10;
salspia_params.delta = 1e-8;
salspia_params.eps_saf = 1e-30;

script_dir = fileparts(mfilename('fullpath'));
repo_dir   = fileparts(script_dir);
addpath(fullfile(repo_dir, 'algorithms'));
addpath(fullfile(repo_dir, 'utils'));

out_dir = fullfile(repo_dir, 'results', 'geometric_error_tables');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% ===== Table A: representative geometric residuals ======================
fprintf('\n============================================================\n');
fprintf('Representative geometric residuals for SaLSPIA\n');
fprintf('============================================================\n');

rep_rows = {};

% Example 4.3: reindeer contour curve
[rms43, max43, it43, einf43] = run_representative_curve_ex43( ...
    p_deg, tol_Ek, maxiter, salspia_params);
rep_rows(end+1, :) = {'Ex.~4.3', 'curve', it43, einf43, rms43, max43}; %#ok<SAGROW>

% Example 4.6: Maungawhau / Mt Eden volcano surface
[rms46, max46, it46, einf46, rel46] = run_representative_volcano_ex46( ...
    p_deg, tol_Ek, maxiter, salspia_params);
rep_rows(end+1, :) = {'Ex.~4.6', 'surface', it46, einf46, rms46, max46}; %#ok<SAGROW>

T_rep = cell2table(rep_rows, 'VariableNames', ...
    {'Example','Type','IT','E_inf','RMSResidual','MaxPointwiseError'});
writetable(T_rep, fullfile(out_dir, 'representative_geometric_errors.csv'));
write_representative_latex(T_rep, fullfile(out_dir, 'representative_geometric_errors.tex'));

fprintf('\nRepresentative geometric residual table:\n');
disp(T_rep);
fprintf('  Volcano relative RMS = %.6f%%\n', rel46);

%% ===== Table B: held-out errors on missing regions ======================
fprintf('\n============================================================\n');
fprintf('Held-out errors on removed data points\n');
fprintf('============================================================\n');

heldout_rows = {};
heldout_rows = append_missing_case(heldout_rows, 'Ex.~4.1', 1, ...
    p_deg, tol_Ek, maxiter, salspia_params, n_repeats);
heldout_rows = append_missing_case(heldout_rows, 'Ex.~4.2', 2, ...
    p_deg, tol_Ek, maxiter, salspia_params, n_repeats);
heldout_rows = append_missing_case(heldout_rows, 'Ex.~4.7', 7, ...
    p_deg, tol_Ek, maxiter, salspia_params, n_repeats);

T_hold = cell2table(heldout_rows, 'VariableNames', ...
    {'Example','Method','ObservedPoints','RemovedPoints','ControlPoints', ...
     'Rank','IT','CPU_mean','CPU_std','E_inf','HoleRMS','HoleMax'});
writetable(T_hold, fullfile(out_dir, 'missing_heldout_errors.csv'));
write_heldout_latex(T_hold, fullfile(out_dir, 'missing_heldout_errors.tex'));

fprintf('\nHeld-out error table (CPU is mean +/- sample SD over five runs):\n');
disp(T_hold(:, {'Example','Method','CPU_mean','CPU_std','HoleRMS','HoleMax'}));

fprintf('\nSaved compact geometric-error outputs to:\n  %s\n', out_dir);

%% ========================================================================
%  Local functions
%% ========================================================================

function [rms_res, max_res, it, einf] = run_representative_curve_ex43( ...
    p_deg, tol_Ek, maxiter, salspia_params)
    m1 = 1000;
    n1 = 288;
    n  = n1 - 1;

    Q = load(salspia_data('cur_data deer'));
    if m1 ~= size(Q, 1)
        idx = round(linspace(1, size(Q, 1), m1));
        Q = Q(idx, :);
    end

    [t, knots] = make_params_knots(Q, n, p_deg);
    A  = build_collocation(knots, p_deg, t);
    P0 = select_initial_ctrl_pts(Q, n);

    [P, info] = salspia(A, Q, P0, tol_Ek, maxiter, salspia_params);
    stat = curve_residual_stats(A, Q, P);

    rms_res = stat.rms_residual;
    max_res = stat.max_pointwise_error;
    it      = info.iter_count;
    einf    = info.Ek_history(end);

    fprintf('  Ex. 4.3 curve: IT=%d, E=%.2e, RMS=%.3e, Max=%.3e\n', ...
        it, einf, rms_res, max_res);
end

function [rms_res, max_res, it, einf, rel_rms_pct] = run_representative_volcano_ex46( ...
    p_deg, tol_Ek, maxiter, salspia_params)
    n1u = 121;
    n1v = 121;

    Z = readmatrix(salspia_data('maungawhau_hr.csv'));
    [m1u, m1v] = size(Z);
    [Qx, Qy] = ndgrid(linspace(0, 865, m1u), linspace(0, 605, m1v));
    Qz = Z;

    u_par = linspace(0, 1, m1u);
    v_par = linspace(0, 1, m1v);
    knots_u = make_clamped_knots_local(n1u, p_deg);
    knots_v = make_clamped_knots_local(n1v, p_deg);
    Au = build_collocation(knots_u, p_deg, u_par);
    Av = build_collocation(knots_v, p_deg, v_par);

    Q_cell = {Qx, Qy, Qz};
    P0_cell = cell(3, 1);
    for c = 1:3
        P0_cell{c} = select_initial_surf_local(Q_cell{c}, n1u, n1v);
    end

    [P, info] = salspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter, salspia_params);
    stat = surface_residual_stats(Au, Av, Q_cell, P);
    elevation_range = max(Qz(:)) - min(Qz(:));

    rms_res = stat.rms_residual;
    max_res = stat.max_pointwise_error;
    it      = info.iter_count;
    einf    = info.Ek_history(end);
    rel_rms_pct = 100 * rms_res / elevation_range;

    fprintf('  Ex. 4.6 volcano: IT=%d, E=%.2e, RMS=%.3e m, Max=%.3e m, RelRMS=%.6f%%\n', ...
        it, einf, rms_res, max_res, rel_rms_pct);
end

function rows = append_missing_case(rows, ex_label, ex_id, ...
    p_deg, tol_Ek, maxiter, salspia_params, n_repeats)
    [A, Q, P0, knots, t_miss, Q_miss] = build_missing_case(ex_id, p_deg);
    A_miss = build_collocation(knots, p_deg, t_miss);

    methods = {'LSPIA-Lin2018', 'ALSPIA', 'SaLSPIA'};
    solvers = {
        @() lspia_lin2018(A, Q, P0, tol_Ek, maxiter), ...
        @() alspia(A, Q, P0, tol_Ek, maxiter), ...
        @() salspia(A, Q, P0, tol_Ek, maxiter, salspia_params)
    };

    fprintf('  %s: observed=%d, removed=%d, control=%d, rank=%d\n', ...
        ex_label, size(Q, 1), size(Q_miss, 1), size(P0, 1), rank(A));

    for im = 1:numel(methods)
        [P, info] = run_with_avg_cpu(solvers{im}, n_repeats);
        C_miss = A_miss * P;
        point_err = sqrt(sum((C_miss - Q_miss).^2, 2));
        hole_rms = sqrt(mean(point_err.^2));
        hole_max = max(point_err);

        fprintf('    %-14s IT=%d, CPU=%.6f +/- %.6f, E=%.2e, Hole RMS=%.3e, Hole Max=%.3e\n', ...
            methods{im}, info.iter_count, info.cpu_time_mean, ...
            info.cpu_time_std, info.Ek_history(end), hole_rms, hole_max);

        rows(end+1, :) = {ex_label, methods{im}, size(Q, 1), size(Q_miss, 1), ...
            size(P0, 1), rank(A), info.iter_count, info.cpu_time_mean, ...
            info.cpu_time_std, info.Ek_history(end), hole_rms, hole_max}; %#ok<AGROW>
    end
end

function [A, Q, P0, knots, t_miss, Q_miss] = build_missing_case(ex_id, p_deg)
    switch ex_id
        case 1
            n = 100;
            m_full = 330;
            theta = linspace(0, 2*pi, m_full+1)';
            r = 2 + 4*cos(2*theta + pi/4) + cos(3*theta + pi/4);
            Q_full = [r .* cos(theta), r .* sin(theta)];
            hole_centers = [0.20, 0.42, 0.62, 0.82];
            hole_width = 0.02;
        case 2
            n = 200;
            m_full = 520;
            theta = linspace(0, 5*6, m_full+1)';
            Q_full = [10*cos(theta*pi/3), 10*sin(theta*pi/3), theta*pi/3];
            hole_centers = [0.25, 0.50, 0.75];
            hole_width = 0.02;
        case 7
            n = 100;
            Q_full = load(salspia_data('s_loop_curve_data.txt'));
            hole_centers = [0.15, 0.38, 0.62, 0.99];
            hole_width = 0.03;
        otherwise
            error('Unknown missing-data example id: %d', ex_id);
    end

    [t_full, knots] = make_params_knots(Q_full, n, p_deg);
    [keep, Q, Q_miss, t_keep, t_miss] = remove_by_param_intervals_local( ...
        Q_full, t_full, hole_centers, hole_width);
    %#ok<ASGLU>
    A  = build_collocation(knots, p_deg, t_keep);
    P0 = select_initial_ctrl_pts(Q_full, n);
end

function [keep, Q_keep, Q_miss, t_keep, t_miss] = remove_by_param_intervals_local( ...
    Q_full, t_full, hole_centers, hole_half_width)
    keep = true(size(t_full));
    for h = 1:length(hole_centers)
        lo = hole_centers(h) - hole_half_width;
        hi = hole_centers(h) + hole_half_width;
        in_hole = (t_full >= lo) & (t_full <= hi);
        keep(in_hole) = false;
    end
    Q_keep = Q_full(keep, :);
    Q_miss = Q_full(~keep, :);
    t_keep = t_full(keep);
    t_miss = t_full(~keep);
end

function knots = make_clamped_knots_local(n1, p)
    n = n1 - 1;
    n_int = n - p;
    knots = zeros(1, n1 + p + 1);
    knots(1:p+1) = 0;
    knots(end-p:end) = 1;
    if n_int > 0
        internal = linspace(0, 1, n_int + 2);
        knots(p+2:end-p-1) = internal(2:end-1);
    end
end

function P0 = select_initial_surf_local(Qc, n1u, n1v)
    [mu1, mv1] = size(Qc);
    idx_u = round(linspace(1, mu1, n1u));
    idx_v = round(linspace(1, mv1, n1v));
    P0 = Qc(idx_u, idx_v);
end

function write_representative_latex(T, fname)
    fid = fopen(fname, 'w');
    if fid == -1
        warning('Cannot write LaTeX table: %s', fname);
        return;
    end
    fprintf(fid, '\\begin{table}[H]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Representative geometric residuals of SaLSPIA. RMS and Max denote the root-mean-square and maximum pointwise residuals, respectively.}\n');
    fprintf(fid, '\\label{tab:representative_geometric_errors}\n');
    fprintf(fid, '\\begin{tabular}{llrr}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'Example & Type & RMS & Max \\\\\n');
    fprintf(fid, '\\hline\n');
    for i = 1:height(T)
        fprintf(fid, '%s & %s & %.3e & %.3e \\\\\n', ...
            T.Example{i}, T.Type{i}, T.RMSResidual(i), T.MaxPointwiseError(i));
    end
    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end

function write_heldout_latex(T, fname)
    fid = fopen(fname, 'w');
    if fid == -1
        warning('Cannot write LaTeX table: %s', fname);
        return;
    end
    fprintf(fid, '\\begin{table}[H]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Held-out geometric errors on removed data points. CPU times are the mean $\\pm$ sample standard deviation over five runs.}\n');
    fprintf(fid, '\\label{tab:missing_heldout_error}\n');
    fprintf(fid, '\\begin{tabular}{llrrrrr}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'Example & Method & IT & CPU (s), mean $\\pm$ SD & $E_\\infty$ & Hole RMS & Hole Max \\\\\n');
    fprintf(fid, '\\hline\n');
    for i = 1:height(T)
        fprintf(fid, '%s & %s & %d & $%.6f \\pm %.6f$ & %.2e & %.3e & %.3e \\\\\n', ...
            T.Example{i}, T.Method{i}, T.IT(i), T.CPU_mean(i), ...
            T.CPU_std(i), T.E_inf(i), T.HoleRMS(i), T.HoleMax(i));
    end
    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end
