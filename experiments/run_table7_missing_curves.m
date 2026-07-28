%% ========================================================================
%  run_table7_missing_curves.m
%  Rank-deficient (missing data) curve fitting experiments
%  (manuscript Table 7, Figs 5-7).
%
%  Four methods compared: LSPIA, LSPIA-Lin2018, ALSPIA, SaLSPIA.
%  Three examples: 4.1 (blob curve), 4.2 (multi-turn helix), 4.7 (G-loop).
%
%  All three missing-data figures share one template (the helix/Figure 6
%  style): missing points shown as red circles, each region tagged with an
%  orange number on the main panel, numbered zoom insets, black fitting
%  curve.
%  ========================================================================
clear; clc; close all;

p_deg   = 3;
tol_Ek  = 1e-6;
maxiter = 10000;
n_repeats = 5;                 % CPU time is averaged over repeated runs

salspia_params.c     = 1e-4;
salspia_params.M     = 10;
salspia_params.delta = 1e-8;
salspia_params.eps_saf = 1e-30;

script_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(script_dir);
out_dir = fullfile(repo_dir, 'results', 'missing_data');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
geom_rows = {};

%% ===== Example 4.1: Blob curve with 4 holes =============================
fprintf('\n============================================================\n');
fprintf('Example 4.1: Blob curve, missing data (4 holes)\n');
fprintf('============================================================\n');

n1 = 100;
m_full1 = 330;
theta1 = linspace(0, 2*pi, m_full1+1)';
r1 = 2 + 4*cos(2*theta1 + pi/4) + cos(3*theta1 + pi/4);
Q_full1 = [r1 .* cos(theta1), r1 .* sin(theta1)];

[t_full1, knots1] = make_params_knots(Q_full1, n1, p_deg);

hole_centers1 = [0.20, 0.42, 0.62, 0.82];
hole_width1   = 0.02;

[keep1, Q1, Q1_miss, t1] = remove_by_param_intervals(...
    Q_full1, t_full1, hole_centers1, hole_width1);

A1 = build_collocation(knots1, p_deg, t1);
P0_1 = select_initial_ctrl_pts(Q_full1, n1);
rk1 = rank(A1);

fprintf('  Data points: %d (removed %d),  Control points: %d,  rank = %d\n', ...
    size(Q1,1), sum(~keep1), n1+1, rk1);

[info_ls1, info_li1, info_al1, info_sl1, P_ls1, P_li1, P_al1, P_sl1] = ...
    run_four_methods(A1, Q1, P0_1, tol_Ek, maxiter, salspia_params, n_repeats);
geom_rows = append_curve_geom_rows(geom_rows, 'Ex4.1', size(Q1,1), n1+1, rk1, ...
    A1, Q1, {'LSPIA','LSPIA-Lin2018','ALSPIA','SaLSPIA'}, ...
    {P_ls1, P_li1, P_al1, P_sl1}, ...
    {info_ls1, info_li1, info_al1, info_sl1}, maxiter);

plot_missing_data_2d(Q1, Q1_miss, P_sl1, knots1, p_deg, ...
    t_full1, hole_centers1, hole_width1);

%% ===== Example 4.2: Helix curve with 3 holes ============================
fprintf('\n============================================================\n');
fprintf('Example 4.2: Helix curve, missing data (3 holes)\n');
fprintf('============================================================\n');

n2 = 200;
m_full2 = 520;

% Multi-turn helix: one full turn needs theta=6, so 5 turns => theta up to 30
theta2 = linspace(0, 5*6, m_full2+1)';
Q_full2 = [10*cos(theta2*pi/3), 10*sin(theta2*pi/3), theta2*pi/3];

[t_full2, knots2] = make_params_knots(Q_full2, n2, p_deg);

hole_centers2 = [0.25, 0.50, 0.75];
hole_width2   = 0.02;

[keep2, Q2, Q2_miss, t2] = remove_by_param_intervals(...
    Q_full2, t_full2, hole_centers2, hole_width2);

A2 = build_collocation(knots2, p_deg, t2);
P0_2 = select_initial_ctrl_pts(Q_full2, n2);
rk2 = rank(A2);

fprintf('  Data points: %d (removed %d),  Control points: %d,  rank = %d\n', ...
    size(Q2,1), sum(~keep2), n2+1, rk2);

[info_ls2, info_li2, info_al2, info_sl2, P_ls2, P_li2, P_al2, P_sl2] = ...
    run_four_methods(A2, Q2, P0_2, tol_Ek, maxiter, salspia_params, n_repeats);
geom_rows = append_curve_geom_rows(geom_rows, 'Ex4.2', size(Q2,1), n2+1, rk2, ...
    A2, Q2, {'LSPIA','LSPIA-Lin2018','ALSPIA','SaLSPIA'}, ...
    {P_ls2, P_li2, P_al2, P_sl2}, ...
    {info_ls2, info_li2, info_al2, info_sl2}, maxiter);

plot_missing_data_3d(Q2, Q2_miss, P_sl2, knots2, p_deg, ...
    t_full2, hole_centers2, hole_width2);

%% ===== Example 4.7: G-loop scattered curve with 4 holes =================
fprintf('\n============================================================\n');
fprintf('Example 4.7: G-loop scattered curve, missing data (4 holes)\n');
fprintf('============================================================\n');

Q_full3 = load(salspia_data('s_loop_curve_data.txt'));
m_full3 = size(Q_full3, 1) - 1;

fprintf('  Loaded %d data points (closed curve)\n', size(Q_full3, 1));

n3 = 100;

[t_full3, knots3] = make_params_knots(Q_full3, n3, p_deg);

hole_centers3 = [0.15, 0.38, 0.62, 0.99];
hole_width3   = 0.03;

[keep3, Q3, Q3_miss, t3] = remove_by_param_intervals(...
    Q_full3, t_full3, hole_centers3, hole_width3);

A3 = build_collocation(knots3, p_deg, t3);
P0_3 = select_initial_ctrl_pts(Q_full3, n3);
rk3 = rank(A3);

fprintf('  Data points: %d (removed %d),  Control points: %d,  rank = %d\n', ...
    size(Q3,1), sum(~keep3), n3+1, rk3);

[info_ls3, info_li3, info_al3, info_sl3, P_ls3, P_li3, P_al3, P_sl3] = ...
    run_four_methods(A3, Q3, P0_3, tol_Ek, maxiter, salspia_params, n_repeats);
geom_rows = append_curve_geom_rows(geom_rows, 'Ex4.7', size(Q3,1), n3+1, rk3, ...
    A3, Q3, {'LSPIA','LSPIA-Lin2018','ALSPIA','SaLSPIA'}, ...
    {P_ls3, P_li3, P_al3, P_sl3}, ...
    {info_ls3, info_li3, info_al3, info_sl3}, maxiter);

plot_missing_data_2d_v2(Q3, Q3_miss, P_sl3, knots3, p_deg, ...
    t_full3, hole_centers3, hole_width3);

%% ===== Convergence History Plot =========================================
figure('Name', 'Convergence: rank-deficient cases', ...
    'Position', [50 50 1600 400]);

examples = {info_ls1, info_li1, info_al1, info_sl1, sprintf('Ex 4.1 (rank=%d)',rk1);
            info_ls2, info_li2, info_al2, info_sl2, sprintf('Ex 4.2 (rank=%d)',rk2);
            info_ls3, info_li3, info_al3, info_sl3, sprintf('Ex 4.7 (rank=%d)',rk3)};

for ie = 1:3
    subplot(1,3,ie);
    ils = examples{ie,1};
    ili = examples{ie,2};
    ial = examples{ie,3};
    isl = examples{ie,4};

    semilogy(0:ils.iter_count, ils.Ek_history, 'b-o', 'MarkerSize', 6, ...
        'MarkerIndices', 1:max(1,floor(max(1,ils.iter_count)/15)):ils.iter_count+1);
    hold on;
    semilogy(0:ili.iter_count, ili.Ek_history, 'm-d', 'MarkerSize', 5, ...
        'MarkerIndices', 1:max(1,floor(max(1,ili.iter_count)/15)):ili.iter_count+1);
    semilogy(0:ial.iter_count, ial.Ek_history, 'r-s', 'MarkerSize', 4, ...
        'MarkerIndices', 1:max(1,floor(max(1,ial.iter_count)/15)):ial.iter_count+1);
    semilogy(0:isl.iter_count, isl.Ek_history, '-p', 'Color', [0.85 0.33 0.10], ...
        'MarkerSize', 5, 'LineWidth', 1.4, ...
        'MarkerIndices', 1:max(1,floor(max(1,isl.iter_count)/15)):isl.iter_count+1);
    hold off;

    yline(tol_Ek, '--', 'Color', [0.5 0.5 0.5]);
    xlabel('Iteration'); ylabel('E_k');
    title(examples{ie,5});
    legend('LSPIA','LSPIA-Lin2018','ALSPIA','SaLSPIA','Location','best');
    grid on;
end

%% ===== Summary Table ====================================================
fprintf('\n\n');
fprintf('=========================================================================================================================================================================\n');
fprintf('                         SUMMARY TABLE: Rank-Deficient Curve Fitting (Missing Data)                                       \n');
fprintf('=========================================================================================================================================================================\n');
fprintf('%-5s %-16s %-6s | %-34s | %-34s | %-34s | %-34s\n', ...
    'Ex', '(m+1,n+1)', 'rank', 'LSPIA', 'LSPIA-Lin2018', 'ALSPIA', 'SaLSPIA');
fprintf('%-5s %-16s %-6s | %-34s | %-34s | %-34s | %-34s\n', ...
    '', '', '', '(IT, CPU mean+/-SD, E_inf)', '(IT, CPU mean+/-SD, E_inf)', ...
    '(IT, CPU mean+/-SD, E_inf)', '(IT, CPU mean+/-SD, E_inf)');
fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');

all_results = {
    info_ls1, info_li1, info_al1, info_sl1, size(Q1,1), n1+1, rk1;
    info_ls2, info_li2, info_al2, info_sl2, size(Q2,1), n2+1, rk2;
    info_ls3, info_li3, info_al3, info_sl3, size(Q3,1), n3+1, rk3
};

for ie = 1:3
    ils = all_results{ie,1};
    ili = all_results{ie,2};
    ial = all_results{ie,3};
    isl = all_results{ie,4};
    m1v = all_results{ie,5};
    n1v = all_results{ie,6};
    rkv = all_results{ie,7};

    if ils.iter_count >= maxiter
        s_ls = '#, #, #';
    else
        s_ls = sprintf('%d, %.6f+/-%0.6f, %.2e', ils.iter_count, ...
            ils.cpu_time_mean, ils.cpu_time_std, ils.Ek_history(end));
    end

    if ili.iter_count >= maxiter
        s_li = '#, #, #';
    else
        s_li = sprintf('%d, %.6f+/-%0.6f, %.2e', ili.iter_count, ...
            ili.cpu_time_mean, ili.cpu_time_std, ili.Ek_history(end));
    end

    if ial.iter_count >= maxiter
        s_al = '#, #, #';
    else
        s_al = sprintf('%d, %.6f+/-%0.6f, %.2e', ial.iter_count, ...
            ial.cpu_time_mean, ial.cpu_time_std, ial.Ek_history(end));
    end

    if isl.iter_count >= maxiter
        s_sl = '#, #, #';
    else
        s_sl = sprintf('%d, %.6f+/-%0.6f, %.2e', isl.iter_count, ...
            isl.cpu_time_mean, isl.cpu_time_std, isl.Ek_history(end));
    end

    fprintf('Ex%-3d (%5d,%5d) %4d  | %-34s | %-34s | %-34s | %-34s\n', ...
        ie, m1v, n1v, rkv, s_ls, s_li, s_al, s_sl);
end
fprintf('=========================================================================================================================================================================\n');
write_missing_geometry_table(geom_rows, out_dir);

function rows = append_curve_geom_rows(rows, ex_label, data_count, ctrl_count, rankA, ...
    A, Q, methods, Ps, infos, maxiter)
    for im = 1:numel(methods)
        info = infos{im};
        stat = curve_residual_stats(A, Q, Ps{im});
        converged = info.iter_count < maxiter;
        if converged
            cpu_mean = info.cpu_time_mean;
            cpu_std = info.cpu_time_std;
            e_inf = info.Ek_history(end);
            rms_residual = stat.rms_residual;
            max_pointwise_error = stat.max_pointwise_error;
            fitting_error = stat.sse;
        else
            cpu_mean = NaN;
            cpu_std = NaN;
            e_inf = NaN;
            rms_residual = NaN;
            max_pointwise_error = NaN;
            fitting_error = NaN;
        end
        rows(end+1, :) = {ex_label, data_count, ctrl_count, rankA, ...
            methods{im}, converged, info.iter_count, cpu_mean, cpu_std, ...
            e_inf, rms_residual, max_pointwise_error, fitting_error}; %#ok<AGROW>
    end
end

function write_missing_geometry_table(rows, out_dir)
    T = cell2table(rows, 'VariableNames', ...
        {'Example','DataPoints','ControlPoints','Rank','Method','Converged', ...
         'IT','CPU_mean','CPU_std','E_inf','RMSResidual', ...
         'MaxPointwiseError','FittingError'});
    writetable(T, fullfile(out_dir, 'missing_data_geometry.csv'));

    fid = fopen(fullfile(out_dir, 'missing_data_geometry.tex'), 'w');
    if fid == -1
        warning('Cannot write missing-data geometry LaTeX table.');
        return;
    end
    fprintf(fid, '\\begin{table}[t]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Geometric residual statistics at the observed data points for rank-deficient curve fitting. CPU times are the mean $\\pm$ sample standard deviation over five runs.}\n');
    fprintf(fid, '\\begin{tabular}{llrrrrr}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'Example & Method & IT & CPU (s), mean $\\pm$ SD & $E_\\infty$ & RMS & Max \\\\\n');
    fprintf(fid, '\\hline\n');
    for i = 1:height(T)
        if T.Converged(i)
            it_text = sprintf('%d', T.IT(i));
            cpu_text = sprintf('$%.6f \\pm %.6f$', T.CPU_mean(i), T.CPU_std(i));
            e_text = sprintf('%.2e', T.E_inf(i));
            rms_text = sprintf('%.3e', T.RMSResidual(i));
            max_text = sprintf('%.3e', T.MaxPointwiseError(i));
        else
            it_text = '\#';
            cpu_text = '\#';
            e_text = '\#';
            rms_text = '\#';
            max_text = '\#';
        end
        fprintf(fid, '%s & %s & %s & %s & %s & %s & %s \\\\\n', ...
            T.Example{i}, T.Method{i}, it_text, cpu_text, e_text, ...
            rms_text, max_text);
    end
    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);

    fprintf('\nSaved missing-data geometric residual outputs to:\n  %s\n', out_dir);
end


%% ========================================================================
%  Helper: remove data by parameter intervals
%% ========================================================================
function [keep, Q_keep, Q_miss, t_keep] = remove_by_param_intervals(...
    Q_full, t_full, hole_centers, hole_half_width)

    keep = true(size(t_full));

    for h = 1:length(hole_centers)
        lo = hole_centers(h) - hole_half_width;
        hi = hole_centers(h) + hole_half_width;
        in_hole = (t_full >= lo) & (t_full <= hi);
        keep(in_hole) = false;

        fprintf('    Hole %d: param [%.4f, %.4f], removed %d pts\n', ...
            h, lo, hi, sum(in_hole));
    end

    fprintf('    Total removed: %d pts, remaining: %d pts\n', ...
        sum(~keep), sum(keep));

    Q_keep = Q_full(keep, :);
    Q_miss = Q_full(~keep, :);
    t_keep = t_full(keep);
end

%% ========================================================================
%  Run four methods
%% ========================================================================
function [info_ls, info_li, info_al, info_sl, ...
          P_ls, P_li, P_al, P_sl] = ...
    run_four_methods(A, Q, P0, tol_Ek, maxiter, salspia_params, n_repeats)

    fprintf('  Running LSPIA (classical) ...');
    [P_ls, info_ls] = run_with_avg_cpu(@() lspia(A, Q, P0, tol_Ek, maxiter), n_repeats);
    if info_ls.iter_count >= maxiter
        fprintf(' FAILED (IT > %d)\n', maxiter);
    else
        fprintf(' IT=%d, CPU=%.6f +/- %.6f, E_inf=%.2e\n', ...
            info_ls.iter_count, info_ls.cpu_time_mean, info_ls.cpu_time_std, ...
            info_ls.Ek_history(end));
    end

    fprintf('  Running LSPIA-Lin2018 ...');
    [P_li, info_li] = run_with_avg_cpu(@() lspia_lin2018(A, Q, P0, tol_Ek, maxiter), n_repeats);
    if info_li.iter_count >= maxiter
        fprintf(' FAILED (IT > %d)\n', maxiter);
    else
        fprintf(' IT=%d, CPU=%.6f +/- %.6f, E_inf=%.2e\n', ...
            info_li.iter_count, info_li.cpu_time_mean, info_li.cpu_time_std, ...
            info_li.Ek_history(end));
    end

    fprintf('  Running ALSPIA ...');
    [P_al, info_al] = run_with_avg_cpu(@() alspia(A, Q, P0, tol_Ek, maxiter), n_repeats);
    fprintf(' IT=%d, CPU=%.6f +/- %.6f, E_inf=%.2e\n', ...
        info_al.iter_count, info_al.cpu_time_mean, info_al.cpu_time_std, ...
        info_al.Ek_history(end));

    fprintf('  Running SaLSPIA ...');
    [P_sl, info_sl] = run_with_avg_cpu(@() salspia(A, Q, P0, tol_Ek, maxiter, salspia_params), n_repeats);
    fprintf(' IT=%d, CPU=%.6f +/- %.6f, E_inf=%.2e\n', ...
        info_sl.iter_count, info_sl.cpu_time_mean, info_sl.cpu_time_std, ...
        info_sl.Ek_history(end));
end

%% ========================================================================
%  red missing circles + orange numbered region labels + numbered insets.
%% ========================================================================
function plot_missing_data_2d(Q, Q_missing, P, knots, p_deg, ...
    t_full, hole_centers, hole_hw)

    t_fine = unique([0, linspace(0, 1, 3000), 1]);
    A_fine = build_collocation(knots, p_deg, t_fine);
    C_fine = A_fine * P;

    holes = compute_hole_insets_from_missing(Q_missing, hole_centers, length(hole_centers));
    n_insets = length(holes);

    col_obs  = [0.40 0.55 0.85];   % observed data : muted blue
    col_miss = [0.85 0.20 0.20];   % missing points: red
    col_fit  = [0.10 0.10 0.10];   % fitting curve : black
    col_lbl  = [0.95 0.45 0.10];   % region number : orange

    % label offset scale (outward from the data centroid)
    ctr = mean(Q, 1);
    sc  = 0.06 * max(max(Q(:,1))-min(Q(:,1)), max(Q(:,2))-min(Q(:,2)));

    figure('Position', [100 100 1400 500]);

    main_w  = 0.30;
    inset_w = 0.08;
    gap = (1.0 - 2*main_w - 2*inset_w) / 5;
    main_bot = 0.13;
    main_h   = 0.78;

    x_main_a   = gap;
    x_insets_a = x_main_a + main_w + gap;
    x_main_b   = x_insets_a + inset_w + gap;
    x_insets_b = x_main_b + main_w + gap;
    inset_h = 0.16;

    % --- (a) Data with missing points ---
    ax_a = axes('Position', [x_main_a, main_bot, main_w, main_h]);
    h_obs = plot(Q(:,1), Q(:,2), '.', 'Color', col_obs, 'MarkerSize', 6); hold on;
    h_mis = plot(Q_missing(:,1), Q_missing(:,2), 'o', 'Color', col_miss, ...
        'MarkerSize', 4, 'LineWidth', 0.8);
    for hi = 1:n_insets
        place_region_label_2d(holes(hi), hi, ctr, sc, col_lbl);
    end
    hold off;
    xlabel('(a) The initial data points with missing data');
    axis equal; grid off; set(gca, 'FontSize', 8);
    legend(ax_a, [h_obs, h_mis], {'Observed data','Missing region'}, ...
        'Location', 'best', 'FontSize', 7);

    for hi = 1:n_insets
        hr = holes(hi);
        inset_y = main_bot + main_h - hi * (inset_h + 0.015);
        axes('Position', [x_insets_a, inset_y, inset_w, inset_h]);
        plot(Q(:,1), Q(:,2), '.', 'Color', col_obs, 'MarkerSize', 5); hold on;
        plot(Q_missing(:,1), Q_missing(:,2), 'o', 'Color', col_miss, ...
            'MarkerSize', 3, 'LineWidth', 0.6);
        hold off;
        xlim([hr.cx - hr.rx, hr.cx + hr.rx]);
        ylim([hr.cy - hr.ry, hr.cy + hr.ry]);
        grid off; set(gca, 'FontSize', 7); box on;
        add_inset_badge(hi, col_lbl);
    end

    % --- (b) Fitting curve ---
    ax_b = axes('Position', [x_main_b, main_bot, main_w, main_h]);
    h_mis2 = plot(Q_missing(:,1), Q_missing(:,2), 'o', 'Color', col_miss, ...
        'MarkerSize', 4, 'LineWidth', 0.8); hold on;
    h_fit = plot(C_fine(:,1), C_fine(:,2), '-', 'Color', col_fit, 'LineWidth', 1.2);
    hold off;
    xlabel('(b) The cubic B-spline fitting curve');
    axis equal; grid off;
    xlim(ax_a.XLim); ylim(ax_a.YLim);
    set(gca, 'FontSize', 8);
    legend(ax_b, [h_mis2, h_fit], {'Missing region','Fitting curve'}, ...
        'Location', 'best', 'FontSize', 7);

    for hi = 1:n_insets
        hr = holes(hi);
        inset_y = main_bot + main_h - hi * (inset_h + 0.015);
        axes('Position', [x_insets_b, inset_y, inset_w, inset_h]);
        plot(Q_missing(:,1), Q_missing(:,2), 'o', 'Color', col_miss, ...
            'MarkerSize', 3, 'LineWidth', 0.6); hold on;
        plot(C_fine(:,1), C_fine(:,2), '-', 'Color', col_fit, 'LineWidth', 1.2);
        hold off;
        xlim([hr.cx - hr.rx, hr.cx + hr.rx]);
        ylim([hr.cy - hr.ry, hr.cy + hr.ry]);
        grid off; set(gca, 'FontSize', 7); box on;
        add_inset_badge(hi, col_lbl);
    end
end

%% ========================================================================
%  red missing circles + orange numbered region labels + numbered insets.
%% ========================================================================
function plot_missing_data_2d_v2(Q, Q_missing, P, knots, p_deg, ...
    t_full, hole_centers, hole_hw)

    t_fine = unique([0, linspace(0, 1, 3000), 1]);
    A_fine = build_collocation(knots, p_deg, t_fine);
    C_fine = A_fine * P;

    holes = compute_hole_insets_tight(Q_missing, hole_centers, length(hole_centers));
    n_insets = length(holes);

    col_obs  = [0.40 0.55 0.85];
    col_miss = [0.85 0.20 0.20];
    col_fit  = [0.10 0.10 0.10];
    col_lbl  = [0.95 0.45 0.10];

    ctr = mean(Q, 1);
    sc  = 0.06 * max(max(Q(:,1))-min(Q(:,1)), max(Q(:,2))-min(Q(:,2)));

    figure('Position', [100 50 1500 600]);

    main_w  = 0.28;
    inset_w = 0.12;
    gap = (1.0 - 2*main_w - 2*inset_w) / 5;
    main_bot = 0.10;
    main_h   = 0.82;

    x_main_a   = gap;
    x_insets_a = x_main_a + main_w + gap;
    x_main_b   = x_insets_a + inset_w + gap;
    x_insets_b = x_main_b + main_w + gap;
    inset_h  = 0.19;
    inset_gap = 0.018;

    % --- (a) Data with missing points ---
    ax_a = axes('Position', [x_main_a, main_bot, main_w, main_h]);
    h_obs = plot(Q(:,1), Q(:,2), '.', 'Color', col_obs, 'MarkerSize', 10); hold on;
    h_mis = plot(Q_missing(:,1), Q_missing(:,2), 'o', 'Color', col_miss, ...
        'MarkerSize', 7, 'LineWidth', 1.0);
    for hi = 1:n_insets
        place_region_label_2d(holes(hi), hi, ctr, sc, col_lbl);
    end
    hold off;
    xlabel('(a) The initial data points with missing data');
    axis equal; grid off; set(gca, 'FontSize', 9);
    legend(ax_a, [h_obs, h_mis], {'Observed data','Missing region'}, ...
        'Location', 'best', 'FontSize', 8);

    for hi = 1:n_insets
        hr = holes(hi);
        inset_y = main_bot + main_h - hi * (inset_h + inset_gap);
        axes('Position', [x_insets_a, inset_y, inset_w, inset_h]);
        plot(Q(:,1), Q(:,2), '.', 'Color', col_obs, 'MarkerSize', 8); hold on;
        plot(Q_missing(:,1), Q_missing(:,2), 'o', 'Color', col_miss, ...
            'MarkerSize', 6, 'LineWidth', 0.8);
        hold off;
        xlim([hr.cx - hr.rx, hr.cx + hr.rx]);
        ylim([hr.cy - hr.ry, hr.cy + hr.ry]);
        grid off; set(gca, 'FontSize', 7); box on;
        add_inset_badge(hi, col_lbl);
    end

    % --- (b) Fitting curve ---
    ax_b = axes('Position', [x_main_b, main_bot, main_w, main_h]);
    h_mis2 = plot(Q_missing(:,1), Q_missing(:,2), 'o', 'Color', col_miss, ...
        'MarkerSize', 7, 'LineWidth', 1.0); hold on;
    h_fit = plot(C_fine(:,1), C_fine(:,2), '-', 'Color', col_fit, 'LineWidth', 1.5);
    hold off;
    xlabel('(b) The cubic B-spline fitting curve');
    axis equal; grid off;
    xlim(ax_a.XLim); ylim(ax_a.YLim);
    set(gca, 'FontSize', 9);
    legend(ax_b, [h_mis2, h_fit], {'Missing region','Fitting curve'}, ...
        'Location', 'best', 'FontSize', 8);

    for hi = 1:n_insets
        hr = holes(hi);
        inset_y = main_bot + main_h - hi * (inset_h + inset_gap);
        axes('Position', [x_insets_b, inset_y, inset_w, inset_h]);
        plot(Q_missing(:,1), Q_missing(:,2), 'o', 'Color', col_miss, ...
            'MarkerSize', 6, 'LineWidth', 0.8); hold on;
        plot(C_fine(:,1), C_fine(:,2), '-', 'Color', col_fit, 'LineWidth', 1.5);
        hold off;
        xlim([hr.cx - hr.rx, hr.cx + hr.rx]);
        ylim([hr.cy - hr.ry, hr.cy + hr.ry]);
        grid off; set(gca, 'FontSize', 7); box on;
        add_inset_badge(hi, col_lbl);
    end
end

%% ========================================================================
%  (Figure 6) 3D helix: distinct missing colour + numbered hole labels
%% ========================================================================
function plot_missing_data_3d(Q, Q_missing, P, knots, p_deg, ...
    t_full, hole_centers, hole_hw)

    t_fine = unique([0, linspace(0, 1, 3000), 1]);
    A_fine = build_collocation(knots, p_deg, t_fine);
    C_fine = A_fine * P;

    holes = compute_helix_hole_insets_by_z(Q_missing, length(hole_centers));
    n_insets = length(holes);
    cents = helix_missing_centroids(Q_missing, length(hole_centers));  % [cx cy cz] per hole

    col_obs  = [0.40 0.55 0.85];
    col_miss = [0.85 0.20 0.20];   % missing points stand out in red for 3D
    col_fit  = [0.10 0.10 0.10];
    col_lbl  = [0.95 0.45 0.10];

    figure('Position', [100 100 1400 500]);

    main_w  = 0.30;
    inset_w = 0.10;
    gap = (1.0 - 2*main_w - 2*inset_w) / 5;
    main_bot = 0.13;
    main_h   = 0.78;

    x_main_a   = gap;
    x_insets_a = x_main_a + main_w + gap;
    x_main_b   = x_insets_a + inset_w + gap;
    x_insets_b = x_main_b + main_w + gap;
    inset_h = 0.22;

    % --- (a) Data with missing points ---
    ax_a = axes('Position', [x_main_a, main_bot, main_w, main_h]);
    h_obs = plot3(Q(:,1), Q(:,2), Q(:,3), '.', 'Color', col_obs, 'MarkerSize', 4); hold on;
    h_mis = plot3(Q_missing(:,1), Q_missing(:,2), Q_missing(:,3), ...
        'o', 'Color', col_miss, 'MarkerSize', 5, 'LineWidth', 0.9);
    for hi = 1:size(cents,1)
        text(cents(hi,1)*1.25, cents(hi,2)*1.25, cents(hi,3), sprintf('%d', hi), ...
            'Color', col_lbl, 'FontSize', 10, 'FontWeight', 'bold', ...
            'BackgroundColor', 'w', 'Margin', 1);
    end
    hold off;
    xlabel('(a) The initial data points with missing data');
    grid off; view(25, 20); set(gca, 'FontSize', 8);
    legend(ax_a, [h_obs, h_mis], {'Observed data','Missing region'}, ...
        'Location', 'best', 'FontSize', 7);

    for hi = 1:n_insets
        hr = holes(hi);
        inset_y = main_bot + main_h - hi * (inset_h + 0.015);
        axes('Position', [x_insets_a, inset_y, inset_w, inset_h]);
        plot3(Q(:,1), Q(:,2), Q(:,3), '.', 'Color', col_obs, 'MarkerSize', 3); hold on;
        plot3(Q_missing(:,1), Q_missing(:,2), Q_missing(:,3), ...
            'o', 'Color', col_miss, 'MarkerSize', 4, 'LineWidth', 0.7);
        hold off;
        xlim([-13, 13]); ylim([-13, 13]);
        zlim([hr.cz - hr.rz, hr.cz + hr.rz]);
        view(25, 20); grid off; set(gca, 'FontSize', 7); box on;
        add_inset_badge(hi, col_lbl);
    end

    % --- (b) Fitting curve ---
    ax_b = axes('Position', [x_main_b, main_bot, main_w, main_h]);
    h_mis2 = plot3(Q_missing(:,1), Q_missing(:,2), Q_missing(:,3), ...
        'o', 'Color', col_miss, 'MarkerSize', 5, 'LineWidth', 0.9); hold on;
    h_fit = plot3(C_fine(:,1), C_fine(:,2), C_fine(:,3), '-', 'Color', col_fit, 'LineWidth', 1.2);
    hold off;
    xlabel('(b) The cubic B-spline fitting curve');
    grid off; view(25, 20);
    xlim(ax_a.XLim); ylim(ax_a.YLim); zlim(ax_a.ZLim);
    set(gca, 'FontSize', 8);
    legend(ax_b, [h_mis2, h_fit], {'Missing region','Fitting curve'}, ...
        'Location', 'best', 'FontSize', 7);

    for hi = 1:n_insets
        hr = holes(hi);
        inset_y = main_bot + main_h - hi * (inset_h + 0.015);
        axes('Position', [x_insets_b, inset_y, inset_w, inset_h]);
        plot3(Q_missing(:,1), Q_missing(:,2), Q_missing(:,3), ...
            'o', 'Color', col_miss, 'MarkerSize', 4, 'LineWidth', 0.7); hold on;
        plot3(C_fine(:,1), C_fine(:,2), C_fine(:,3), '-', 'Color', col_fit, 'LineWidth', 1.2);
        hold off;
        xlim([-13, 13]); ylim([-13, 13]);
        zlim([hr.cz - hr.rz, hr.cz + hr.rz]);
        view(25, 20); grid off; set(gca, 'FontSize', 7); box on;
        add_inset_badge(hi, col_lbl);
    end
end

%% ========================================================================
%  Compute inset regions: by Z coordinate - for multi-turn helix (Ex 4.2)
%% ========================================================================
function holes = compute_helix_hole_insets_by_z(Q_missing, n_holes)

    if isempty(Q_missing)
        holes = struct('cz',{},'rz',{});
        return;
    end

    z = Q_missing(:, 3);
    z_sorted = sort(z);
    N = length(z_sorted);

    % Find large gaps in z to split into clusters
    dz = diff(z_sorted);
    n_gaps = n_holes - 1;
    if n_gaps > 0 && N > n_holes
        [~, gap_idx] = sort(dz, 'descend');
        split_pts = sort(gap_idx(1:min(n_gaps, length(gap_idx))));
        cluster_starts = [1; split_pts+1];
        cluster_ends   = [split_pts; N];
    else
        cluster_starts = 1;
        cluster_ends   = N;
    end

    holes = struct('cz',{},'rz',{});
    for c = 1:length(cluster_starts)
        ci = cluster_starts(c):cluster_ends(c);
        Z  = z_sorted(ci);
        cz = mean(Z);
        % Half one helix turn height = pi, add small margin
        rz = max(abs(Z - cz)) + pi + 0.5;
        rz = max(rz, pi + 0.5);
        holes(end+1).cz = cz;
        holes(end).rz   = rz;
    end
end

%% ========================================================================
%  Compute inset regions: ORIGINAL loose padding (for Ex 4.1 blob)
%% ========================================================================
function holes = compute_hole_insets_from_missing(Q_missing, hole_centers, n_expected)

    if isempty(Q_missing) || n_expected == 0
        holes = struct('cx',{},'cy',{},'rx',{},'ry',{});
        return;
    end

    N = size(Q_missing, 1);
    dists = sqrt(sum(diff(Q_missing(:,1:2)).^2, 2));

    n_gaps = n_expected - 1;
    if n_gaps > 0 && N > n_expected
        [~, gap_idx] = sort(dists, 'descend');
        split_points = sort(gap_idx(1:min(n_gaps, length(gap_idx))));
        cluster_starts = [1; split_points + 1];
        cluster_ends   = [split_points; N];
    else
        cluster_starts = 1;
        cluster_ends = N;
    end

    holes = struct('cx',{},'cy',{},'rx',{},'ry',{});

    for c = 1:length(cluster_starts)
        ci = cluster_starts(c):cluster_ends(c);
        if length(ci) < 1, continue; end

        X = Q_missing(ci, 1);
        Y = Q_missing(ci, 2);

        cx = mean(X);
        cy = mean(Y);
        rx = max(abs(X - cx)) * 2.0 + 0.2;
        ry = max(abs(Y - cy)) * 2.0 + 0.2;
        rx = max(rx, 0.4);
        ry = max(ry, 0.4);

        holes(end+1).cx = cx;
        holes(end).cy = cy;
        holes(end).rx = rx;
        holes(end).ry = ry;
    end
end

%% ========================================================================
%  Compute inset regions: TIGHT padding (for Ex 4.7 G-loop)
%% ========================================================================
function holes = compute_hole_insets_tight(Q_missing, hole_centers, n_expected)

    if isempty(Q_missing) || n_expected == 0
        holes = struct('cx',{},'cy',{},'rx',{},'ry',{});
        return;
    end

    N = size(Q_missing, 1);
    dists = sqrt(sum(diff(Q_missing(:,1:2)).^2, 2));

    n_gaps = n_expected - 1;
    if n_gaps > 0 && N > n_expected
        [~, gap_idx] = sort(dists, 'descend');
        split_points = sort(gap_idx(1:min(n_gaps, length(gap_idx))));
        cluster_starts = [1; split_points + 1];
        cluster_ends   = [split_points; N];
    else
        cluster_starts = 1;
        cluster_ends = N;
    end

    holes = struct('cx',{},'cy',{},'rx',{},'ry',{});

    for c = 1:length(cluster_starts)
        ci = cluster_starts(c):cluster_ends(c);
        if length(ci) < 1, continue; end

        X = Q_missing(ci, 1);
        Y = Q_missing(ci, 2);

        cx = mean(X);
        cy = mean(Y);
        % Tighter: 1.5x spread + 0.05 padding, min 0.08
        rx = max(abs(X - cx)) * 1.5 + 0.05;
        ry = max(abs(Y - cy)) * 1.5 + 0.05;
        rx = max(rx, 0.08);
        ry = max(ry, 0.08);

        holes(end+1).cx = cx;
        holes(end).cy = cy;
        holes(end).rx = rx;
        holes(end).ry = ry;
    end
end

%% ========================================================================
%  Helper: orange numbered label for a 2D missing region, placed just
%  outside the region (radially outward from the data centroid)
%% ========================================================================
function place_region_label_2d(hr, idx, ctr, sc, col)
    d = [hr.cx, hr.cy] - ctr;
    if norm(d) < eps, d = [1 0]; end
    d = d / norm(d);
    lx = hr.cx + sc * d(1);
    ly = hr.cy + sc * d(2);
    text(lx, ly, sprintf('%d', idx), 'Color', col, ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'BackgroundColor', 'w', 'Margin', 1);
end

%% ========================================================================
%  Helper: small numbered badge in the top-left corner of an inset
%% ========================================================================
function add_inset_badge(idx, col)
    text(0.06, 0.90, sprintf('%d', idx), 'Units', 'normalized', ...
        'Color', col, 'FontSize', 8, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'Margin', 1, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end

%% ========================================================================
%  Helper: centroid [cx cy cz] of each helix hole's missing points
%% ========================================================================
function cents = helix_missing_centroids(Q_missing, n_holes)
    cents = zeros(0,3);
    if isempty(Q_missing), return; end

    z = Q_missing(:,3);
    [z_sorted, order] = sort(z);
    Qs = Q_missing(order, :);
    N = numel(z_sorted);

    dz = diff(z_sorted);
    n_gaps = n_holes - 1;
    if n_gaps > 0 && N > n_holes
        [~, gap_idx] = sort(dz, 'descend');
        split_pts = sort(gap_idx(1:min(n_gaps, numel(gap_idx))));
        cluster_starts = [1; split_pts+1];
        cluster_ends   = [split_pts; N];
    else
        cluster_starts = 1;
        cluster_ends   = N;
    end

    for c = 1:numel(cluster_starts)
        ci = cluster_starts(c):cluster_ends(c);
        cents(end+1, :) = mean(Qs(ci, :), 1); %#ok<AGROW>
    end
end
