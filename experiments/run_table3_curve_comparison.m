%% ========================================================================
%  run_table3_curve_comparison.m
%  Curve fitting comparison (paper Table 3, Examples 4.1-4.3):
%      LSPIA  vs  ALSPIA  vs  MLSPIA  vs  NmLSPIA  vs  SaLSPIA
%
%  SaLSPIA = salspia.m  (spectral-adaptive LSPIA, the proposed method).
%  Run setup_paths.m once before this script.
%  ========================================================================
clear; clc; close all;

%% ===== Settings =========================================================
p_deg   = 3;
tol_Ek  = 1e-6;
maxiter = 10000;
n_repeats = 5;                 % CPU time is averaged over repeated runs

salspia_params.c     = 1e-4;   % GLL sufficient-decrease (sigma)
salspia_params.M     = 10;     % GLL memory length
salspia_params.delta = 1e-8;   % truncation threshold (eq 8)
salspia_params.eps_saf = 1e-30; % safeguard tolerance (Algorithm 1)

% (example_id, m+1, n+1)
configs = {
    1, 15001, 3001;
    2, 10001, 2001;
    3,  1000,  288;
};

%% ===== Run ==============================================================
n_cfg   = size(configs, 1);
results = cell(n_cfg, 1);

for ic = 1:n_cfg
    ex_id = configs{ic,1};  m1 = configs{ic,2};  n1 = configs{ic,3};
    m = m1-1;  n = n1-1;

    fprintf('\n=== Example %d:  m+1=%d,  n+1=%d ===\n', ex_id, m1, n1);

    Q           = generate_curve_data(ex_id, m1);
    [t, knots]  = make_params_knots(Q, n, p_deg);
    A           = build_collocation(knots, p_deg, t);
    P0          = select_initial_ctrl_pts(Q, n);

    eigvals = sort(real(eig(A'*A)));
    eigvals = eigvals(eigvals > 1e-14);
    kappa   = eigvals(end)/eigvals(1);
    fprintf('  kappa = %.2e\n', kappa);

    fprintf('  LSPIA      ... ');
    [P_ls, info_ls] = run_with_avg_cpu(@() lspia(A, Q, P0, tol_Ek, maxiter), n_repeats);
    fprintf('IT=%d  CPU=%.6f +/- %.6fs  E=%.2e\n', ...
        info_ls.iter_count, info_ls.cpu_time_mean, info_ls.cpu_time_std, info_ls.Ek_history(end));

    fprintf('  ALSPIA     ... ');
    [P_al, info_al] = run_with_avg_cpu(@() alspia(A, Q, P0, tol_Ek, maxiter), n_repeats);
    fprintf('IT=%d  CPU=%.6f +/- %.6fs  E=%.2e\n', ...
        info_al.iter_count, info_al.cpu_time_mean, info_al.cpu_time_std, info_al.Ek_history(end));

    fprintf('  MLSPIA     ... ');
    [P_ml, info_ml] = run_with_avg_cpu(@() mlspia(A, Q, P0, tol_Ek, maxiter), n_repeats);
    fprintf('IT=%d  CPU=%.6f +/- %.6fs  E=%.2e  (omega=%.4f, nu=%.4f)\n', ...
        info_ml.iter_count, info_ml.cpu_time_mean, info_ml.cpu_time_std, ...
        info_ml.Ek_history(end), info_ml.omega, info_ml.nu);

    fprintf('  NmLSPIA    ... ');
    [P_nm, info_nm] = run_with_avg_cpu(@() nmlspia(A, Q, P0, tol_Ek, maxiter), n_repeats);
    fprintf('IT=%d  CPU=%.6f +/- %.6fs  E=%.2e  (zeta=%.4f, eta=%.4f)\n', ...
        info_nm.iter_count, info_nm.cpu_time_mean, info_nm.cpu_time_std, ...
        info_nm.Ek_history(end), info_nm.zeta, info_nm.eta);

    fprintf('  SaLSPIA    ... ');
    [P_sl, info_sl] = run_with_avg_cpu(@() salspia(A, Q, P0, tol_Ek, maxiter, salspia_params), n_repeats);
    fprintf('IT=%d  CPU=%.6f +/- %.6fs  E=%.2e\n', ...
        info_sl.iter_count, info_sl.cpu_time_mean, info_sl.cpu_time_std, info_sl.Ek_history(end));

    stat_ls = curve_residual_stats(A, Q, P_ls);
    stat_al = curve_residual_stats(A, Q, P_al);
    stat_ml = curve_residual_stats(A, Q, P_ml);
    stat_nm = curve_residual_stats(A, Q, P_nm);
    stat_sl = curve_residual_stats(A, Q, P_sl);

    res.ex_id=ex_id; res.m1=m1; res.n1=n1; res.kappa=kappa;
    res.info_ls=info_ls; res.info_al=info_al; res.info_ml=info_ml;
    res.info_nm=info_nm; res.info_sl=info_sl;
    res.stat_ls=stat_ls; res.stat_al=stat_al; res.stat_ml=stat_ml;
    res.stat_nm=stat_nm; res.stat_sl=stat_sl;
    results{ic} = res;

    %% Plot
    figure('Name',sprintf('Curve Ex%d (m=%d,n=%d)',ex_id,m,n),'Position',[100 100 1200 450]);
    subplot(1,2,1);
    mk_ls = 1:max(1,floor(info_ls.iter_count/15)):info_ls.iter_count+1;
    mk_al = 1:max(1,floor(info_al.iter_count/12)):info_al.iter_count+1;
    mk_ml = 1:max(1,floor(info_ml.iter_count/12)):info_ml.iter_count+1;
    mk_nm = 1:max(1,floor(info_nm.iter_count/12)):info_nm.iter_count+1;
    mk_sl = 1:max(1,floor(info_sl.iter_count/12)):info_sl.iter_count+1;

    col_nm = [0.47 0.67 0.19];   % NmLSPIA green

    semilogy(0:info_ls.iter_count, info_ls.Ek_history,'b-o','MarkerSize',3,'DisplayName','LSPIA','MarkerIndices',mk_ls);
    hold on;
    semilogy(0:info_al.iter_count, info_al.Ek_history,'r-s','MarkerSize',4,'DisplayName','ALSPIA','MarkerIndices',mk_al);
    semilogy(0:info_ml.iter_count, info_ml.Ek_history,'m-d','MarkerSize',4,'DisplayName','MLSPIA','MarkerIndices',mk_ml);
    semilogy(0:info_nm.iter_count, info_nm.Ek_history,'-^','Color',col_nm,'MarkerSize',4,'DisplayName','NmLSPIA','MarkerIndices',mk_nm);
    semilogy(0:info_sl.iter_count, info_sl.Ek_history,'-p','Color',[0.85 0.33 0.10],'MarkerSize',5,'LineWidth',1.4,'DisplayName','SaLSPIA','MarkerIndices',mk_sl);
    hold off;
    yline(tol_Ek,'--','Color',[.5 .5 .5],'HandleVisibility','off');
    xlabel('Iteration k'); ylabel('E_k');
    title(sprintf('Ex %d: m+1=%d, n+1=%d  (\\kappa=%.1e)',ex_id,m1,n1,kappa));
    legend('Location','best'); grid on;

    subplot(1,2,2);
    semilogy(0:info_ls.iter_count, info_ls.err_history,'b-o','MarkerSize',3,'DisplayName','LSPIA','MarkerIndices',mk_ls);
    hold on;
    semilogy(0:info_al.iter_count, info_al.err_history,'r-s','MarkerSize',4,'DisplayName','ALSPIA','MarkerIndices',mk_al);
    semilogy(0:info_ml.iter_count, info_ml.err_history,'m-d','MarkerSize',4,'DisplayName','MLSPIA','MarkerIndices',mk_ml);
    semilogy(0:info_nm.iter_count, info_nm.err_history,'-^','Color',col_nm,'MarkerSize',4,'DisplayName','NmLSPIA','MarkerIndices',mk_nm);
    semilogy(0:info_sl.iter_count, info_sl.err_history,'-p','Color',[0.85 0.33 0.10],'MarkerSize',5,'LineWidth',1.4,'DisplayName','SaLSPIA','MarkerIndices',mk_sl);
    hold off;
    xlabel('Iteration k'); ylabel('||AP^k-Q||_F^2');
    title('Fitting error'); legend('Location','best'); grid on;
    drawnow;
end

%% ===== Summary Table ====================================================
labels = {'LSPIA','ALSPIA','MLSPIA','NmLSPIA','SaLSPIA'};
fields = {'info_ls','info_al','info_ml','info_nm','info_sl'};
fprintf('\n');
fprintf('=========================================================================================================================================================\n');
fprintf('  SUMMARY TABLE (Curve Fitting)   --   E_inf = E_k at termination\n');
fprintf('---------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%-6s %-14s %9s', 'Ex', '(m+1,n+1)', 'kappa');
for j = 1:numel(labels), fprintf(' | %-30s', labels{j}); end
fprintf('\n%-6s %-14s %9s', '', '', '');
for j = 1:numel(labels), fprintf(' | %6s %14s %8s', 'IT','CPU mean+/-SD','E_inf'); end
fprintf('\n---------------------------------------------------------------------------------------------------------------------------------------------------------\n');
for ic = 1:n_cfg
    r = results{ic};
    fprintf('Ex %d   (%5d,%5d) %9.1e', r.ex_id, r.m1, r.n1, r.kappa);
    for j = 1:numel(fields)
        info = r.(fields{j});
        fprintf(' | %6d %6.4f+/-%6.4f %8.2e', info.iter_count, ...
            info.cpu_time_mean, info.cpu_time_std, info.Ek_history(end));
    end
    fprintf('\n');
end

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'curve_comparison');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
write_curve_geometry_tables(results, out_dir);


%% ===== Data Generator ===================================================
function Q = generate_curve_data(ex_id, m1)
    switch ex_id
        case 1
            theta = linspace(0,2*pi,m1)';
            r = 2 + 4*cos(2*theta+pi/4) + cos(3*theta+pi/4);
            Q = [r.*cos(theta), r.*sin(theta)];
        case 2
            theta = linspace(0,2*pi,m1)';
            Q = [10*cos(theta*pi/3), 10*sin(theta*pi/3), theta*pi/3];
        case 3
            fname = salspia_data('cur_data deer');
            if exist(fname,'file')
                Q = load(fname);
                if m1 ~= size(Q,1)
                    idx = round(linspace(1,size(Q,1),m1));
                    Q   = Q(idx,:);
                end
            else
                warning('Reindeer data not found; using synthetic substitute.');
                theta = linspace(0,2*pi,m1)';
                r = 1+0.3*sin(7*theta)+0.15*cos(13*theta);
                Q = [r.*cos(theta), r.*sin(theta)];
            end
        otherwise
            error('Unknown example id: %d', ex_id);
    end
end

function write_curve_geometry_tables(results, out_dir)
    labels = {'LSPIA','ALSPIA','MLSPIA','NmLSPIA','SaLSPIA'};
    info_fields = {'info_ls','info_al','info_ml','info_nm','info_sl'};
    stat_fields = {'stat_ls','stat_al','stat_ml','stat_nm','stat_sl'};
    rows = {};
    for ic = 1:numel(results)
        r = results{ic};
        for im = 1:numel(labels)
            info = r.(info_fields{im});
            stat = r.(stat_fields{im});
            rows(end+1, :) = {r.ex_id, r.m1, r.n1, r.kappa, labels{im}, ...
                info.iter_count, info.cpu_time_mean, info.cpu_time_std, ...
                info.Ek_history(end), ...
                stat.rms_residual, stat.max_pointwise_error, stat.sse}; %#ok<AGROW>
        end
    end

    T = cell2table(rows, 'VariableNames', ...
        {'Example','DataPoints','ControlPoints','Kappa','Method','IT', ...
         'CPU_mean','CPU_std','E_inf','RMSResidual','MaxPointwiseError', ...
         'FittingError'});
    writetable(T, fullfile(out_dir, 'curve_comparison_geometry.csv'));
    write_curve_cpu_latex(T, fullfile(out_dir, 'table3_curve_cpu.tex'));

    fid = fopen(fullfile(out_dir, 'curve_comparison_geometry.tex'), 'w');
    if fid == -1
        warning('Cannot write curve geometry LaTeX table.');
        return;
    end
    fprintf(fid, '\\begin{table}[t]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Geometric residual statistics for curve fitting in Examples~4.1--4.3.}\n');
    fprintf(fid, '\\begin{tabular}{llrrrr}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'Example & Method & IT & $E_\\infty$ & RMS & Max \\\\\n');
    fprintf(fid, '\\hline\n');
    for i = 1:height(T)
        fprintf(fid, 'Ex.~%d & %s & %d & %.2e & %.3e & %.3e \\\\\n', ...
            T.Example(i), T.Method{i}, T.IT(i), T.E_inf(i), ...
            T.RMSResidual(i), T.MaxPointwiseError(i));
    end
    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);

    fprintf('\nSaved curve geometric residual outputs to:\n  %s\n', out_dir);
end

function write_curve_cpu_latex(T, fname)
    fid = fopen(fname, 'w');
    if fid == -1
        warning('Cannot write curve CPU LaTeX table.');
        return;
    end

    labels = {'LSPIA','ALSPIA','MLSPIA','NmLSPIA','SaLSPIA'};
    examples = [1, 2, 3];

    fprintf(fid, '\\begin{table}[t]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Numerical results for curve fitting. CPU times are the mean $\\pm$ sample standard deviation over five runs.}\n');
    fprintf(fid, '\\label{tab:curve_results_cpu_sd}\n');
    fprintf(fid, '\\begin{tabular}{llccc}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'Method & Metric & $(15001,3001)$ & $(10001,2001)$ & $(1000,288)$ \\\\\n');
    fprintf(fid, '\\hline\n');

    for im = 1:numel(labels)
        rows = cell(1, numel(examples));
        for ie = 1:numel(examples)
            idx = T.Example == examples(ie) & strcmp(T.Method, labels{im});
            rows{ie} = T(idx, :);
        end

        fprintf(fid, '%s & $E_\\infty$ & %.2e & %.2e & %.2e \\\\\n', ...
            labels{im}, rows{1}.E_inf, rows{2}.E_inf, rows{3}.E_inf);
        fprintf(fid, ' & IT & %d & %d & %d \\\\\n', ...
            rows{1}.IT, rows{2}.IT, rows{3}.IT);
        fprintf(fid, ' & CPU (s) & $%.4f \\pm %.4f$ & $%.4f \\pm %.4f$ & $%.4f \\pm %.4f$ \\\\\n', ...
            rows{1}.CPU_mean, rows{1}.CPU_std, rows{2}.CPU_mean, ...
            rows{2}.CPU_std, rows{3}.CPU_mean, rows{3}.CPU_std);
    end

    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end
