%% ========================================================================
%  run_table5_surface_comparison.m
%  Surface fitting comparison (paper Table 5, Examples 4.4-4.5):
%      LSPIA  vs  ALSPIA  vs  MLSPIA  vs  NmLSPIA  vs  SaLSPIA
%
%  Examples 4.4 (Dini) and 4.5 (peaks) are analytic (no data file).
%  SaLSPIA = salspia_surf.m (the proposed method).
%  Run setup_paths.m once before this script.
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

% (ex_id, m1, m2, n1, n2)
configs = {
    4,  81,  81, 31, 31;
    5, 121, 121, 41, 41;
};

n_cfg   = size(configs,1);
results = cell(n_cfg,1);

for ic = 1:n_cfg
    ex_id=configs{ic,1}; m1u=configs{ic,2}; m1v=configs{ic,3};
    n1u=configs{ic,4};   n1v=configs{ic,5};

    fprintf('\n=== Example %d: data=%dx%d, control=%dx%d ===\n',...
        ex_id, m1u, m1v, n1u, n1v);

    [Qx,Qy,Qz,u_raw,v_raw] = generate_surface_data(ex_id, m1u, m1v);
    u_par = (u_raw-u_raw(1))/(u_raw(end)-u_raw(1));
    v_par = (v_raw-v_raw(1))/(v_raw(end)-v_raw(1));

    knots_u = make_clamped_knots(n1u, p_deg);
    knots_v = make_clamped_knots(n1v, p_deg);

    Au = build_collocation(knots_u, p_deg, u_par);
    Av = build_collocation(knots_v, p_deg, v_par);

    Q_cell  = {Qx, Qy, Qz};
    P0_cell = cell(3,1);
    for c=1:3, P0_cell{c} = select_initial_surf(Q_cell{c}, n1u, n1v); end

    fprintf('  LSPIA      ... ');
    [P_ls, info_ls] = run_with_avg_cpu(@() lspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter), n_repeats);
    fprintf('IT=%d  CPU=%.6f +/- %.6fs  E=%.2e\n', ...
        info_ls.iter_count, info_ls.cpu_time_mean, info_ls.cpu_time_std, info_ls.Ek_history(end));

    fprintf('  ALSPIA     ... ');
    [P_al, info_al] = run_with_avg_cpu(@() alspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter), n_repeats);
    fprintf('IT=%d  CPU=%.6f +/- %.6fs  E=%.2e\n', ...
        info_al.iter_count, info_al.cpu_time_mean, info_al.cpu_time_std, info_al.Ek_history(end));

    fprintf('  MLSPIA     ... ');
    [P_ml, info_ml] = run_with_avg_cpu(@() mlspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter), n_repeats);
    fprintf('IT=%d  CPU=%.6f +/- %.6fs  E=%.2e  (omega=%.4f, nu=%.6f)\n',...
        info_ml.iter_count, info_ml.cpu_time_mean, info_ml.cpu_time_std, ...
        info_ml.Ek_history(end), info_ml.omega, info_ml.nu);

    fprintf('  NmLSPIA    ... ');
    [P_nm, info_nm] = run_with_avg_cpu(@() nmlspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter), n_repeats);
    fprintf('IT=%d  CPU=%.6f +/- %.6fs  E=%.2e  (zeta=%.6f, eta=%.6f)\n',...
        info_nm.iter_count, info_nm.cpu_time_mean, info_nm.cpu_time_std, ...
        info_nm.Ek_history(end), info_nm.zeta, info_nm.eta);

    fprintf('  SaLSPIA    ... ');
    [P_sl, info_sl] = run_with_avg_cpu(@() salspia_surf(Au, Av, Q_cell, P0_cell, tol_Ek, maxiter, salspia_params), n_repeats);
    fprintf('IT=%d  CPU=%.6f +/- %.6fs  E=%.2e\n', ...
        info_sl.iter_count, info_sl.cpu_time_mean, info_sl.cpu_time_std, info_sl.Ek_history(end));

    stat_ls = surface_residual_stats(Au, Av, Q_cell, P_ls);
    stat_al = surface_residual_stats(Au, Av, Q_cell, P_al);
    stat_ml = surface_residual_stats(Au, Av, Q_cell, P_ml);
    stat_nm = surface_residual_stats(Au, Av, Q_cell, P_nm);
    stat_sl = surface_residual_stats(Au, Av, Q_cell, P_sl);

    res.ex_id=ex_id; res.m1u=m1u; res.m1v=m1v; res.n1u=n1u; res.n1v=n1v;
    res.info_ls=info_ls; res.info_al=info_al; res.info_ml=info_ml; res.info_nm=info_nm; res.info_sl=info_sl;
    res.stat_ls=stat_ls; res.stat_al=stat_al; res.stat_ml=stat_ml; res.stat_nm=stat_nm; res.stat_sl=stat_sl;
    results{ic} = res;

    %% Plot
    figure('Name',sprintf('Surface Ex%d (%d,%d,%d,%d)',ex_id,m1u,m1v,n1u,n1v),...
           'Position',[100 100 700 500]);
    mk_ls = 1:max(1,floor(info_ls.iter_count/20)):info_ls.iter_count+1;
    mk_al = 1:max(1,floor(info_al.iter_count/15)):info_al.iter_count+1;
    mk_ml = 1:max(1,floor(info_ml.iter_count/15)):info_ml.iter_count+1;
    mk_nm = 1:max(1,floor(info_nm.iter_count/15)):info_nm.iter_count+1;
    mk_sl = 1:max(1,floor(info_sl.iter_count/15)):info_sl.iter_count+1;

    semilogy(0:info_ls.iter_count, info_ls.Ek_history,'b-o','MarkerSize',3,'DisplayName','LSPIA','MarkerIndices',mk_ls);
    hold on;
    semilogy(0:info_al.iter_count, info_al.Ek_history,'r-s','MarkerSize',4,'DisplayName','ALSPIA','MarkerIndices',mk_al);
    semilogy(0:info_ml.iter_count, info_ml.Ek_history,'m-d','MarkerSize',4,'DisplayName','MLSPIA','MarkerIndices',mk_ml);
    semilogy(0:info_nm.iter_count, info_nm.Ek_history,'g-^','MarkerSize',4,'DisplayName','NmLSPIA','MarkerIndices',mk_nm);
    semilogy(0:info_sl.iter_count, info_sl.Ek_history,'-p','Color',[0.85 0.33 0.10],'MarkerSize',5,'LineWidth',1.4,'DisplayName','SaLSPIA','MarkerIndices',mk_sl);
    hold off;
    yline(tol_Ek,'--','Color',[.5 .5 .5],'HandleVisibility','off');
    xlabel('Iteration k'); ylabel('E_k');
    title(sprintf('Ex %d: data=%d\\times%d, control=%d\\times%d',ex_id,m1u,m1v,n1u,n1v));
    legend('Location','best'); grid on;
    drawnow;
end

%% Summary
labels = {'LSPIA','ALSPIA','MLSPIA','NmLSPIA','SaLSPIA'};
fields = {'info_ls','info_al','info_ml','info_nm','info_sl'};
fprintf('\n');
fprintf('=================================================================================================================================================================\n');
fprintf('  SUMMARY TABLE (Surface Fitting)\n');
fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%-6s %-18s', 'Ex', '(dataU,dataV,ctrlU,ctrlV)');
for j = 1:numel(labels), fprintf(' | %-30s', labels{j}); end
fprintf('\n%-6s %-18s', '', '');
for j = 1:numel(labels), fprintf(' | %6s %14s %8s', 'IT','CPU mean+/-SD','E_inf'); end
fprintf('\n-----------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
for ic = 1:n_cfg
    r = results{ic};
    fprintf('Ex %d  (%3d,%3d,%3d,%3d)', r.ex_id, r.m1u, r.m1v, r.n1u, r.n1v);
    for j = 1:numel(fields)
        info = r.(fields{j});
        fprintf(' | %6d %6.4f+/-%6.4f %8.2e', info.iter_count, ...
            info.cpu_time_mean, info.cpu_time_std, info.Ek_history(end));
    end
    fprintf('\n');
end

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'surface_comparison');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
write_surface_geometry_tables(results, out_dir);

%% Helpers
function [Qx,Qy,Qz,u_par,v_par] = generate_surface_data(ex_id, m1u, m1v)
    switch ex_id
        case 4
            th1=linspace(0,4*pi,m1u); th2=linspace(0.1,pi-0.1,m1v);
            [T1,T2]=meshgrid(th1,th2); T1=T1'; T2=T2';
            Qx=cos(T1).*sin(T2);
            Qy=sin(T1).*sin(T2);
            Qz=cos(T2)+log(tan(T2/2))+T1;
            u_par=th1(:); v_par=th2(:);
        case 5
            th1=linspace(-3,3,m1u); th2=linspace(-4,4,m1v);
            [T1,T2]=meshgrid(th1,th2); T1=T1'; T2=T2';
            F=3*(1-T1).^2.*exp(-T1.^2-(T2+1).^2) ...
              -10*(T1/5-T1.^3-T2.^5).*exp(-T1.^2-T2.^2) ...
              -1/3*exp(-(T1+1).^2-T2.^2);
            Qx=T1; Qy=T2; Qz=F;
            u_par=th1(:); v_par=th2(:);
        otherwise
            error('Unknown surface example: %d',ex_id);
    end
end

function knots = make_clamped_knots(n1, p)
    n=n1-1; n_int=n-p;
    knots=zeros(1,n1+p+1);
    knots(1:p+1)=0; knots(end-p:end)=1;
    if n_int>0
        internal=linspace(0,1,n_int+2);
        knots(p+2:end-p-1)=internal(2:end-1);
    end
end

function P0 = select_initial_surf(Qc, n1u, n1v)
    [mu1,mv1]=size(Qc);
    idx_u=round(linspace(1,mu1,n1u));
    idx_v=round(linspace(1,mv1,n1v));
    P0=Qc(idx_u,idx_v);
end

function write_surface_geometry_tables(results, out_dir)
    labels = {'LSPIA','ALSPIA','MLSPIA','NmLSPIA','SaLSPIA'};
    info_fields = {'info_ls','info_al','info_ml','info_nm','info_sl'};
    stat_fields = {'stat_ls','stat_al','stat_ml','stat_nm','stat_sl'};
    rows = {};
    for ic = 1:numel(results)
        r = results{ic};
        for im = 1:numel(labels)
            info = r.(info_fields{im});
            stat = r.(stat_fields{im});
            rows(end+1, :) = {r.ex_id, r.m1u, r.m1v, r.n1u, r.n1v, ...
                labels{im}, info.iter_count, info.cpu_time_mean, ...
                info.cpu_time_std, info.Ek_history(end), stat.rms_residual, ...
                stat.max_pointwise_error, stat.sse}; %#ok<AGROW>
        end
    end

    T = cell2table(rows, 'VariableNames', ...
        {'Example','DataU','DataV','CtrlU','CtrlV','Method','IT', ...
         'CPU_mean','CPU_std','E_inf','RMSResidual','MaxPointwiseError', ...
         'FittingError'});
    writetable(T, fullfile(out_dir, 'surface_comparison_geometry.csv'));
    write_surface_cpu_latex(T, fullfile(out_dir, 'table5_surface_ex44_ex45_cpu.tex'));

    fid = fopen(fullfile(out_dir, 'surface_comparison_geometry.tex'), 'w');
    if fid == -1
        warning('Cannot write surface geometry LaTeX table.');
        return;
    end
    fprintf(fid, '\\begin{table}[t]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Geometric residual statistics for surface fitting in Examples~4.4--4.5.}\n');
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
end

function write_surface_cpu_latex(T, fname)
    fid = fopen(fname, 'w');
    if fid == -1
        warning('Cannot write surface CPU LaTeX table.');
        return;
    end

    labels = {'LSPIA','ALSPIA','MLSPIA','NmLSPIA','SaLSPIA'};
    examples = [4, 5];

    fprintf(fid, '\\begin{table}[t]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Numerical results for surface fitting in Examples~4.4--4.5. CPU times are the mean $\\pm$ sample standard deviation over five runs.}\n');
    fprintf(fid, '\\label{tab:surface_results_cpu_sd_part1}\n');
    fprintf(fid, '\\begin{tabular}{llcc}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'Method & Metric & $(81,81,31,31)$ & $(121,121,41,41)$ \\\\\n');
    fprintf(fid, '\\hline\n');

    for im = 1:numel(labels)
        rows = cell(1, numel(examples));
        for ie = 1:numel(examples)
            idx = T.Example == examples(ie) & strcmp(T.Method, labels{im});
            rows{ie} = T(idx, :);
        end

        fprintf(fid, '%s & $E_\\infty$ & %.2e & %.2e \\\\\n', ...
            labels{im}, rows{1}.E_inf, rows{2}.E_inf);
        fprintf(fid, ' & IT & %d & %d \\\\\n', rows{1}.IT, rows{2}.IT);
        fprintf(fid, ' & CPU (s) & $%.4f \\pm %.4f$ & $%.4f \\pm %.4f$ \\\\\n', ...
            rows{1}.CPU_mean, rows{1}.CPU_std, ...
            rows{2}.CPU_mean, rows{2}.CPU_std);
    end

    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end
