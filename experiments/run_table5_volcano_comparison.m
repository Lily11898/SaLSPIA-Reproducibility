%% ========================================================================
%  run_table5_volcano_comparison.m
%  Volcano comparison for Table 5 and iteration-progression Figure 4
%  (Example 4.6, Maungawhau LiDAR).
%  Comparison: LSPIA vs ALSPIA vs MLSPIA vs NmLSPIA vs SaLSPIA; the iteration
%  snapshots in the figure visualise SaLSPIA (= salspia_surf).
%
%  Data file: maungawhau_hr.csv, resolved from ../data
%    860 x 600 elevation samples for the Maungawhau LiDAR example.
%
%  Run setup_paths.m once before this script.
%  ========================================================================
clear; clc; close all;

%% ===== User Settings ====================================================
script_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(script_dir);

p_deg   = 3;
tol_Ek  = 1e-6;
maxiter = 10000;
n_repeats = 5;                 % CPU time is averaged over repeated runs

salspia_params.c     = 1e-4;
salspia_params.M     = 10;
salspia_params.delta = 1e-8;
salspia_params.eps_saf = 1e-30;

n1u = 121;
n1v = 121;

data_file = 'maungawhau_hr.csv';

out_dir = fullfile(repo_dir, 'results', 'volcano_surface_fit');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ===== Load Data ========================================================
fprintf('Loading data from: %s\n', data_file);
Z = readmatrix(salspia_data(data_file));
[m1u, m1v] = size(Z);

fprintf('  Grid: (m1+1, m2+1) = (%d, %d)\n', m1u, m1v);
fprintf('  Total vertices: %d\n', m1u * m1v);
fprintf('  Control points: (n1+1, n2+1) = (%d, %d)\n', n1u, n1v);

[Qx, Qy] = ndgrid(linspace(0, 865, m1u), linspace(0, 605, m1v));
Qz = Z;

%% ===== Setup B-spline ===================================================
u_par = linspace(0, 1, m1u);
v_par = linspace(0, 1, m1v);

knots_u = make_clamped_knots(n1u, p_deg);
knots_v = make_clamped_knots(n1v, p_deg);

fprintf('Building collocation matrices ... ');
tic;
Au = build_collocation(knots_u, p_deg, u_par);
Av = build_collocation(knots_v, p_deg, v_par);
t_build = toc;
fprintf('done (%.3f s)\n', t_build);

Q_cell  = {Qx, Qy, Qz};
P0_cell = cell(3, 1);
for c = 1:3
    P0_cell{c} = select_initial_surf(Q_cell{c}, n1u, n1v);
end

%% ===== Run Five Algorithms ===============================================
fprintf('\n--- Running LSPIA ---\n');
[~, info_ls] = run_with_avg_cpu(@() lspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter), n_repeats);
fprintf('  IT = %d,  CPU = %.6f +/- %.6f s,  E_inf = %.2e\n', ...
    info_ls.iter_count, info_ls.cpu_time_mean, info_ls.cpu_time_std, info_ls.Ek_history(end));

fprintf('--- Running ALSPIA ---\n');
[~, info_al] = run_with_avg_cpu(@() alspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter), n_repeats);
fprintf('  IT = %d,  CPU = %.6f +/- %.6f s,  E_inf = %.2e\n', ...
    info_al.iter_count, info_al.cpu_time_mean, info_al.cpu_time_std, info_al.Ek_history(end));

fprintf('--- Running MLSPIA ---\n');
[~, info_ml] = run_with_avg_cpu(@() mlspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter), n_repeats);
fprintf('  IT = %d,  CPU = %.6f +/- %.6f s,  E_inf = %.2e\n', ...
    info_ml.iter_count, info_ml.cpu_time_mean, info_ml.cpu_time_std, info_ml.Ek_history(end));
fprintf('  (omega = %.6f,  nu = %.6f)\n', info_ml.omega, info_ml.nu);

fprintf('--- Running NmLSPIA ---\n');
[~, info_nm] = run_with_avg_cpu(@() nmlspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter), n_repeats);
fprintf('  IT = %d,  CPU = %.6f +/- %.6f s,  E_inf = %.2e\n', ...
    info_nm.iter_count, info_nm.cpu_time_mean, info_nm.cpu_time_std, info_nm.Ek_history(end));
fprintf('  (zeta = %.6f,  eta = %.6f)\n', info_nm.zeta, info_nm.eta);

fprintf('--- Running SaLSPIA ---\n');
[P_sl, info_sl] = run_with_avg_cpu(@() salspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter, salspia_params), n_repeats);
fprintf('  IT = %d,  CPU = %.6f +/- %.6f s,  E_inf = %.2e\n', ...
    info_sl.iter_count, info_sl.cpu_time_mean, info_sl.cpu_time_std, info_sl.Ek_history(end));

%% ===== Summary Table ====================================================
fprintf('\n');
fprintf('==============================================================================\n');
fprintf('  SUMMARY:  grid = %dx%d,  ctrl = %dx%d\n', m1u, m1v, n1u, n1v);
fprintf('====================================================================================================\n');
fprintf('  %-15s | %8s  %12s  %12s  %12s\n', ...
    'Algorithm', 'IT', 'CPU mean(s)', 'CPU SD(s)', 'E_inf');
fprintf('  %-15s | %8s  %12s  %12s  %12s\n', ...
    '---------------', '--------', '------------', '------------', '------------');
fprintf('  %-15s | %8d  %12.6f  %12.6f  %12.2e\n', 'LSPIA', info_ls.iter_count, ...
    info_ls.cpu_time_mean, info_ls.cpu_time_std, info_ls.Ek_history(end));
fprintf('  %-15s | %8d  %12.6f  %12.6f  %12.2e\n', 'ALSPIA', info_al.iter_count, ...
    info_al.cpu_time_mean, info_al.cpu_time_std, info_al.Ek_history(end));
fprintf('  %-15s | %8d  %12.6f  %12.6f  %12.2e\n', 'MLSPIA', info_ml.iter_count, ...
    info_ml.cpu_time_mean, info_ml.cpu_time_std, info_ml.Ek_history(end));
fprintf('  %-15s | %8d  %12.6f  %12.6f  %12.2e\n', 'NmLSPIA', info_nm.iter_count, ...
    info_nm.cpu_time_mean, info_nm.cpu_time_std, info_nm.Ek_history(end));
fprintf('  %-15s | %8d  %12.6f  %12.6f  %12.2e\n', 'SaLSPIA', info_sl.iter_count, ...
    info_sl.cpu_time_mean, info_sl.cpu_time_std, info_sl.Ek_history(end));

method = {'LSPIA'; 'ALSPIA'; 'MLSPIA'; 'NmLSPIA'; 'SaLSPIA'};
info_all = {info_ls; info_al; info_ml; info_nm; info_sl};
IT = cellfun(@(s) s.iter_count, info_all);
CPU_mean = cellfun(@(s) s.cpu_time_mean, info_all);
CPU_std = cellfun(@(s) s.cpu_time_std, info_all);
E_inf = cellfun(@(s) s.Ek_history(end), info_all);
cpu_summary = table(method, IT, CPU_mean, CPU_std, E_inf, ...
    'VariableNames', {'Method','IT','CPU_mean','CPU_std','E_inf'});
writetable(cpu_summary, fullfile(out_dir, 'volcano_cpu_summary.csv'));
write_volcano_cpu_latex(cpu_summary, fullfile(out_dir, 'table5_surface_ex46_cpu.tex'));
save(fullfile(out_dir, 'volcano_cpu_results.mat'), ...
    'cpu_summary', 'info_ls', 'info_al', 'info_ml', 'info_nm', 'info_sl');



%% ===== Local residual statistics for Table 6 ===========================
Sx_sl = Au * P_sl{1} * Av';
Sy_sl = Au * P_sl{2} * Av';
Sz_sl = Au * P_sl{3} * Av';
Rerr = sqrt((Sx_sl - Qx).^2 + (Sy_sl - Qy).^2 + (Sz_sl - Qz).^2);
rms_residual = sqrt(mean(Rerr(:).^2));
elevation_range = max(Qz(:)) - min(Qz(:));
relative_rms_pct = 100 * rms_residual / elevation_range;
max_residual = max(Rerr(:));

fprintf('\nLocal residual statistics for SaLSPIA (raw pointwise residuals):\n');
fprintf('  RMS residual    = %.6e m\n', rms_residual);
fprintf('  Elevation range = %.6e m\n', elevation_range);
fprintf('  Relative RMS    = %.6f%%\n', relative_rms_pct);
fprintf('  Max residual    = %.6e m\n', max_residual);

residual_summary = table(rms_residual, elevation_range, relative_rms_pct, ...
    max_residual, 'VariableNames', {'RMSResidual', 'ElevationRange', ...
    'RelativeRMSPercent', 'MaxResidual'});
writetable(residual_summary, fullfile(out_dir, 'volcano_residual_summary.csv'));

%% ========================================================================
%  PLOT:  3-panel iteration progression
%  Layout:  [a. iter 0]  [b. iter mid]  [c. limit]
%
%   * Surface coloured by elevation (hypsometric tints), shared colour
%     range across all 3 panels so they are directly comparable.
%   * Subsampled control polygon overlaid as small dark dots - shows how
%     the control net is reshaped during the iteration.
%
%  Visualised algorithm: SaLSPIA (manuscript Figure 4).
% =========================================================================

%% ----- Configurable -----------------------------------------------------
N_PANEL     = 3;                 % 3 or 4 panels
iter_choice = 'auto';            % 'auto' or 'fixed'
iter_fixed  = [0, 7, 13];        % used if iter_choice = 'fixed'
                                 %   For the paper figure, the automatic choice
                                 %   gives start / mid / converged snapshots.
                                 %   Switch iter_choice to 'auto' if the total
                                 %   IT changes on a different dataset.
view_angle  = [-35, 50];         % [azimuth, elevation]
                                 %   alternatives: [-45, 55] or [-25, 60]

% Vertical exaggeration (use 1.0 for true scale, smaller = more exaggerated)
z_aspect = 0.5;                  % daspect([1 1 z_aspect])

% Algorithm visualised in the iteration-progression figure: SaLSPIA.
algo_name = 'SaLSPIA';
algo_run  = @(maxit) salspia_surf(Au, Av, Q_cell, P0_cell, ...
                                  1e-30, maxit, salspia_params);
total_it  = info_sl.iter_count;

%% ----- Decide iteration list -------------------------------------------
if strcmpi(iter_choice, 'auto')
    if N_PANEL == 3
        iter_list = [0, max(1, round(total_it/2)), total_it];
    else
        iter_list = [0, max(1, round(total_it/3)), ...
                        max(2, round(2*total_it/3)), total_it];
    end
else
    iter_list = iter_fixed;
end
assert(numel(iter_list) == N_PANEL, ...
    'iter_list length (%d) must match N_PANEL (%d).', numel(iter_list), N_PANEL);
fprintf('\nComputing iteration snapshots for %s: [%s ] ...\n', ...
        algo_name, sprintf(' %d', iter_list));

%% ----- Compute control points + E_inf at each snapshot -----------------
P_snap  = cell(N_PANEL, 1);
Ek_snap = zeros(N_PANEL, 1);
for s = 1:N_PANEL
    if iter_list(s) == 0
        P_snap{s} = P0_cell;
        % Initial E_inf at iter 0  (max pointwise residual norm)
        Sx0 = Au * P0_cell{1} * Av';
        Sy0 = Au * P0_cell{2} * Av';
        Sz0 = Au * P0_cell{3} * Av';
        Ek_snap(s) = max(sqrt((Sx0(:) - Qx(:)).^2 + ...
                              (Sy0(:) - Qy(:)).^2 + ...
                              (Sz0(:) - Qz(:)).^2));
    else
        [P_snap{s}, info_s] = algo_run(iter_list(s));
        Ek_snap(s) = info_s.Ek_history(end);
    end
end

%% ----- High-resolution evaluation grid ---------------------------------
N_eval = 200;
u_eval = linspace(0, 1, N_eval);
v_eval = linspace(0, 1, N_eval);
Au_ev  = build_collocation(knots_u, p_deg, u_eval);
Av_ev  = build_collocation(knots_v, p_deg, v_eval);

%% ----- Shared colour range (from data) ---------------------------------
z_lo = min(Qz(:));
z_hi = max(Qz(:));
dem_cmap = make_dem_colormap(256);

%% ----- Show ALL control points -----------------------------------------
ctrl_idx_u = 1 : n1u;
ctrl_idx_v = 1 : n1v;

%% ----- Build the figure ------------------------------------------------
fig_w = 480 * N_PANEL + 80;        % auto-size width based on N_PANEL
fig_h = 620;                       % a bit taller to make room for (a)(b)(c) labels
fig = figure('Name', sprintf('%s iteration progression', algo_name), ...
             'Position', [40 120 fig_w fig_h], 'Color', 'w');

ax_handles = gobjects(N_PANEL, 1);
for s = 1:N_PANEL
    ax_handles(s) = subplot(1, N_PANEL, s);

    % Evaluate fitted surface at this snapshot (dense grid)
    Sx = Au_ev * P_snap{s}{1} * Av_ev';
    Sy = Au_ev * P_snap{s}{2} * Av_ev';
    Sz = Au_ev * P_snap{s}{3} * Av_ev';

    % Surface coloured by elevation
    surf(Sx, Sy, Sz, Sz, ...
         'EdgeColor',         'none', ...
         'FaceLighting',      'gouraud', ...
         'AmbientStrength',   0.55, ...
         'DiffuseStrength',   0.85, ...
         'SpecularStrength',  0.10, ...
         'SpecularExponent',  8);
    hold on;

    % Control polygon nodes, lifted slightly
    Pcx = P_snap{s}{1}(ctrl_idx_u, ctrl_idx_v);
    Pcy = P_snap{s}{2}(ctrl_idx_u, ctrl_idx_v);
    Pcz = P_snap{s}{3}(ctrl_idx_u, ctrl_idx_v);
    z_lift = 0.005 * (z_hi - z_lo);
    plot3(Pcx(:), Pcy(:), Pcz(:) + z_lift, ...
          '.', 'Color', [0.23 0.04 0.47], 'MarkerSize', 3);

    hold off;

    colormap(ax_handles(s), dem_cmap);
    caxis(ax_handles(s), [z_lo, z_hi]);   %#ok<CAXIS>
    camlight('headlight');
    material([0.55 0.85 0.10 8 1]);
    axis tight; axis off;
    daspect([1 1 z_aspect]);
    view(view_angle(1), view_angle(2));
end

%% ----- Panel labels: (a) k=0, (b) k=7, (c) k=15 ------------------------
panel_letters = {'(a)', '(b)', '(c)', '(d)'};
for s = 1:N_PANEL
    pos_s = get(ax_handles(s), 'Position');
    annotation('textbox', ...
        [pos_s(1), pos_s(2) - 0.02, pos_s(3), 0.05], ...
        'String',              sprintf('%s $k = %d$', panel_letters{s}, iter_list(s)), ...
        'EdgeColor',           'none', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'top', ...
        'Interpreter',         'latex', ...
        'FontSize',            18, ...
        'Margin',              1);
end

%% ----- One shared colourbar on the right side --------------------------
last_pos = get(ax_handles(end), 'Position');     % [x y w h]
cb_h = colorbar(ax_handles(end));
set(cb_h, 'Position', [last_pos(1) + last_pos(3) + 0.012, ...
                       last_pos(2) + 0.10, ...
                       0.010, ...
                       last_pos(4) - 0.20]);
ylabel(cb_h, 'Elevation (m)', 'FontSize', 11);

fprintf('Done. %d-panel iteration figure generated for %s.\n', N_PANEL, algo_name);
fprintf('Saved volcano residual outputs to:\n  %s\n', out_dir);

%% ===== Helper Functions =================================================
function knots = make_clamped_knots(n1, p)
    n          = n1 - 1;
    n_internal = n - p;
    knots      = zeros(1, n1 + p + 1);
    knots(1:p+1)       = 0;
    knots(end-p:end)   = 1;
    if n_internal > 0
        internal = linspace(0, 1, n_internal + 2);
        knots(p+2 : end-p-1) = internal(2:end-1);
    end
end

function P0 = select_initial_surf(Qc, n1u, n1v)
    [mu1, mv1] = size(Qc);
    idx_u = round(linspace(1, mu1, n1u));
    idx_v = round(linspace(1, mv1, n1v));
    P0    = Qc(idx_u, idx_v);
end

function cmap = make_dem_colormap(n)
% MAKE_DEM_COLORMAP   Hypsometric tints for terrain visualisation.
%   Geographic colour ramp: low ground = green, mid = tan / brown,
%   peaks = light grey.  Pure RGB anchors, interpolated with pchip.
    if nargin < 1, n = 256; end
    anchors = [
        0.00  0.20 0.45 0.20;   % dark green  (low / valley)
        0.20  0.45 0.65 0.30;   % medium green
        0.40  0.85 0.80 0.50;   % yellow-tan
        0.60  0.75 0.55 0.30;   % brown
        0.80  0.60 0.40 0.25;   % dark brown
        1.00  0.95 0.95 0.95];  % light grey  (summit)
    t = linspace(0, 1, n)';
    r = interp1(anchors(:,1), anchors(:,2), t, 'pchip');
    g = interp1(anchors(:,1), anchors(:,3), t, 'pchip');
    b = interp1(anchors(:,1), anchors(:,4), t, 'pchip');
    cmap = max(0, min(1, [r, g, b]));
end

function write_volcano_cpu_latex(T, fname)
fid = fopen(fname, 'w');
if fid == -1
    warning('Cannot write volcano CPU LaTeX table.');
    return;
end

fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Numerical results for the volcano surface in Example~4.6. CPU times are the mean $\\pm$ sample standard deviation over five runs.}\n');
fprintf(fid, '\\label{tab:surface_results_cpu_sd_ex46}\n');
fprintf(fid, '\\begin{tabular}{lrrr}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Method & $E_\\infty$ & IT & CPU (s), mean $\\pm$ SD \\\\\n');
fprintf(fid, '\\hline\n');
for i = 1:height(T)
    fprintf(fid, '%s & %.2e & %d & $%.4f \\pm %.4f$ \\\\\n', ...
        T.Method{i}, T.E_inf(i), T.IT(i), T.CPU_mean(i), T.CPU_std(i));
end
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
end
