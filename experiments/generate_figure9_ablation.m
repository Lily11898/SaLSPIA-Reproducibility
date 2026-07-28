%% ========================================================================
%  generate_figure9_ablation.m
%  Convergence histories for the three weight configurations in Example 4.7.
%
%  This figure supports the ablation discussion:
%    1. Long endpoint  : alpha_k^(1) = BB1 endpoint, with GLL safeguard
%    2. Short endpoint : alpha_k^(2) = BB2 endpoint, with GLL safeguard
%    3. SaLSPIA        : interpolated spectral weight, with GLL safeguard
%
%  Output:
%    results/ablation_convergence_ex47/fig_ablation_convergence_ex47_envelope.pdf
%    results/ablation_convergence_ex47/fig_ablation_convergence_ex47_envelope.png
%    results/ablation_convergence_ex47/ablation_convergence_ex47_summary.csv
%
%  The envelope figure plots log10(min_{0<=j<=k} E_j). It is often clearer
%  for the missing-data case because the GLL criterion is nonmonotone.
%  ========================================================================
clear; clc; close all;

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(this_file));
addpath(fullfile(repo_root, 'algorithms'));
addpath(fullfile(repo_root, 'utils'));

p_deg   = 3;
tol_Ek  = 1e-6;
maxiter = 10000;

params.c           = 1e-4;
params.M           = 10;
params.delta       = 1e-8;
params.eps_saf     = 1e-30;
params.p           = 3;
params.e_tol       = 1e-15;
params.max_halving = 100;

out_dir = fullfile(repo_root, 'results', 'ablation_convergence_ex47');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ===== Build Example 4.7 missing-data system ===========================
Q_full = load(salspia_data('s_loop_curve_data.txt'));
n = 100;

[t_full, knots] = make_params_knots(Q_full, n, p_deg);

hole_centers = [0.15, 0.38, 0.62, 0.99];
hole_width   = 0.03;

[keep, Q_obs, ~, t_obs] = remove_by_param_intervals( ...
    Q_full, t_full, hole_centers, hole_width);

A = build_collocation(knots, p_deg, t_obs);
P0 = select_initial_ctrl_pts(Q_full, n);
rk = rank(A);

fprintf('\nExample 4.7 ablation convergence plot\n');
fprintf('  Observed data points: %d\n', size(Q_obs, 1));
fprintf('  Removed data points : %d\n', sum(~keep));
fprintf('  Control points      : %d\n', n + 1);
fprintf('  rank(A)             : %d\n\n', rk);

%% ===== Run the three weight configurations =============================
methods = {
    'Long endpoint $\alpha_k^{(1)}$',  'bb1';
    'Short endpoint $\alpha_k^{(2)}$', 'bb2';
    'SaLSPIA',                         'salspia';
};

infos = cell(size(methods, 1), 1);

for im = 1:size(methods, 1)
    run_params = params;
    run_params.mode = methods{im, 2};

    fprintf('  %-16s ... ', methods{im, 2});
    [~, info] = ablation_lspia(A, Q_obs, P0, tol_Ek, maxiter, run_params);
    infos{im} = info;

    fprintf('IT=%d, halvings=%d, E_final=%.3e, status=%s\n', ...
        info.iter_count, info.total_halving, info.Ek_history(end), info.status);
end

%% ===== Save summary data ===============================================
summary = cell(size(methods, 1), 6);
for im = 1:size(methods, 1)
    info = infos{im};
    summary(im, :) = {methods{im, 1}, methods{im, 2}, info.iter_count, ...
        info.total_halving, info.Ek_history(end), info.status};
end

T = cell2table(summary, 'VariableNames', ...
    {'Method', 'Mode', 'IT', 'Halvings', 'E_final', 'Status'});
writetable(T, fullfile(out_dir, 'ablation_convergence_ex47_summary.csv'));
disp(T);

%% ===== Export manuscript Figure 9 ======================================
plot_envelope(infos, methods, tol_Ek, out_dir);

fprintf('\nSaved outputs to:\n  %s\n', out_dir);

%% ===== Helpers ==========================================================
function [keep, Q_keep, Q_miss, t_keep] = remove_by_param_intervals( ...
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
end

function plot_envelope(infos, methods, tol_Ek, out_dir)
fig = figure('Color', 'w', 'Position', [100 100 620 410]);
ax = axes(fig);
hold(ax, 'on');

line_styles = {'-', '--', '-'};
markers = {'o', 's', '^'};
colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.10 0.60 0.20];

for im = 1:size(methods, 1)
    info = infos{im};
    y = max(info.Ek_history(:), realmin);
    y = cummin(y);
    kvec = 0:info.iter_count;
    ylog = log10(y);
    marker_step = max(2, ceil(numel(kvec) / 8));

    plot(ax, kvec, ylog, ...
        'LineStyle', line_styles{im}, ...
        'Marker', markers{im}, ...
        'Color', colors(im, :), ...
        'LineWidth', 1.9, ...
        'MarkerSize', 5.5, ...
        'MarkerFaceColor', 'w', ...
        'MarkerIndices', 1:marker_step:numel(kvec), ...
        'DisplayName', methods{im, 1});
end

yline(ax, log10(tol_Ek), ':', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.2, 'DisplayName', 'Tolerance');

hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.16;
ax.MinorGridAlpha = 0.08;
ax.FontName = 'Times New Roman';
ax.FontSize = 12;

xmax = max(cellfun(@(s) s.iter_count, infos));
xlim(ax, [0, xmax]);
ylim(ax, [-7, 0.15]);
set(ax, 'YTick', -7:1:0);

xlabel(ax, 'Iteration number $k$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel(ax, '$\log_{10}(\min_{0\leq j\leq k}E_j)$', ...
    'Interpreter', 'latex', 'FontSize', 13);
legend(ax, 'Location', 'southwest', 'Interpreter', 'latex', ...
    'FontSize', 11, 'Box', 'off');

base = fullfile(out_dir, 'fig_ablation_convergence_ex47_envelope');
safe_export(fig, [base '.pdf']);
safe_export(fig, [base '.png']);
end

function safe_export(fig, fname)
try
    if endsWith(fname, '.pdf')
        exportgraphics(fig, fname, 'ContentType', 'vector');
    else
        exportgraphics(fig, fname, 'Resolution', 600);
    end
catch
    saveas(fig, fname);
end
end
