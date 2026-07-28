%% ========================================================================
%  run_table8_ablation.m
%  Formal reproduction entry for manuscript Table 8.
%
%  Compares:
%    - long endpoint  alpha_k^(1);
%    - short endpoint alpha_k^(2);
%    - SaLSPIA interpolated weight;
%  on Curve Ex. 4.2, Curve Ex. 4.3, and missing Curve Ex. 4.7.
%
%  Outputs:
%    results/table8_ablation/table8_ablation.csv
%    results/table8_ablation/table8_ablation.tex
%  ========================================================================
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
repo_dir   = fileparts(script_dir);
addpath(fullfile(repo_dir, 'algorithms'));
addpath(fullfile(repo_dir, 'utils'));

p_deg = 3;
tol_Ek = 1e-6;
maxiter = 10000;

params.c = 1e-4;
params.M = 10;
params.delta = 1e-8;
params.eps_saf = 1e-30;
params.p = 3;
params.e_tol = 1e-15;
params.max_halving = 100;

methods = {
    'Long endpoint alpha^(1)',  'bb1';
    'Short endpoint alpha^(2)', 'bb2';
    'SaLSPIA',                  'salspia'
};

cases = build_cases(p_deg);
rows = {};

fprintf('\n============================================================\n');
fprintf('Ablation of the spectrally adaptive weight (Table 8)\n');
fprintf('============================================================\n');

for ic = 1:numel(cases)
    cs = cases(ic);
    fprintf('\n%s: data=%d, control=%d, rank=%d\n', ...
        cs.label, size(cs.A,1), size(cs.A,2), cs.rankA);
    for im = 1:size(methods,1)
        run_params = params;
        run_params.mode = methods{im,2};
        [~, info] = ablation_lspia(cs.A, cs.Q, cs.P0, ...
            tol_Ek, maxiter, run_params);
        fprintf('  %-28s IT=%d, Halvings=%d, E=%.2e\n', ...
            methods{im,1}, info.iter_count, info.total_halving, ...
            info.Ek_history(end));
        rows(end+1,:) = {cs.label, size(cs.A,1), size(cs.A,2), ...
            cs.rankA, methods{im,1}, methods{im,2}, ...
            info.Ek_history(end), info.iter_count, ...
            info.total_halving}; %#ok<SAGROW>
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'Example','DataPoints','ControlPoints','Rank','Method','Mode', ...
     'E_inf','IT','Halvings'});

out_dir = fullfile(repo_dir, 'results', 'table8_ablation');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
writetable(T, fullfile(out_dir, 'table8_ablation.csv'));
write_table8_latex(T, fullfile(out_dir, 'table8_ablation.tex'));

fprintf('\nSaved Table 8 outputs to:\n  %s\n', out_dir);

function cases = build_cases(p_deg)
cases = struct('label',{},'A',{},'Q',{},'P0',{},'rankA',{});

% Curve Ex. 4.2: analytic helix.
m1 = 10001; n1 = 2001; n = n1 - 1;
theta = linspace(0, 2*pi, m1)';
Q = [10*cos(theta*pi/3), 10*sin(theta*pi/3), theta*pi/3];
[t, knots] = make_params_knots(Q, n, p_deg);
A = build_collocation(knots, p_deg, t);
cases(1) = make_case('Curve Ex.~4.2', A, Q, ...
    select_initial_ctrl_pts(Q, n));

% Curve Ex. 4.3: reindeer contour, resampled to 1000 points.
m1 = 1000; n1 = 288; n = n1 - 1;
Q = load(salspia_data('cur_data deer'));
if size(Q,1) ~= m1
    idx = round(linspace(1, size(Q,1), m1));
    Q = Q(idx,:);
end
[t, knots] = make_params_knots(Q, n, p_deg);
A = build_collocation(knots, p_deg, t);
cases(2) = make_case('Curve Ex.~4.3', A, Q, ...
    select_initial_ctrl_pts(Q, n));

% Curve Ex. 4.7: G-shaped loop with four missing intervals.
n = 100;
Q_full = load(salspia_data('s_loop_curve_data.txt'));
[t_full, knots] = make_params_knots(Q_full, n, p_deg);
keep = true(size(t_full));
centers = [0.15, 0.38, 0.62, 0.99];
half_width = 0.03;
for h = 1:numel(centers)
    keep(t_full >= centers(h)-half_width & ...
         t_full <= centers(h)+half_width) = false;
end
Q = Q_full(keep,:);
A = build_collocation(knots, p_deg, t_full(keep));
cases(3) = make_case('Curve Ex.~4.7', A, Q, ...
    select_initial_ctrl_pts(Q_full, n));
end

function cs = make_case(label, A, Q, P0)
cs.label = label;
cs.A = A;
cs.Q = Q;
cs.P0 = P0;
cs.rankA = rank(A);
end

function write_table8_latex(T, fname)
fid = fopen(fname, 'w');
if fid == -1, error('Cannot write LaTeX table: %s', fname); end
cleanup = onCleanup(@() fclose(fid));

examples = {'Curve Ex.~4.2','Curve Ex.~4.3','Curve Ex.~4.7'};
method_labels = {
    'Long endpoint $\alpha_k^{(1)}$';
    'Short endpoint $\alpha_k^{(2)}$';
    'SaLSPIA'
};
modes = {'bb1','bb2','salspia'};

fprintf(fid, '\\begin{table}[H]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Ablation of the spectrally adaptive weight on three curve examples.}\n');
fprintf(fid, '\\label{tab:abl_weight}\n');
fprintf(fid, '\\begin{tabular}{llccc}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Method & & Curve Ex.~4.2 & Curve Ex.~4.3 & Curve Ex.~4.7 \\\\\n');
fprintf(fid, '& data/control & $(10001,2001)$ & $(1000,288)$ & $(206,101)$ \\\\\n');
fprintf(fid, '& rank & full & full & 90 \\\\\n');
fprintf(fid, '\\hline\n');

for im = 1:numel(modes)
    idx1 = strcmp(T.Example, examples{1}) & strcmp(T.Mode, modes{im});
    idx2 = strcmp(T.Example, examples{2}) & strcmp(T.Mode, modes{im});
    idx3 = strcmp(T.Example, examples{3}) & strcmp(T.Mode, modes{im});
    fprintf(fid, '%s & $E_\\infty$ & %.2e & %.2e & %.2e \\\\\n', ...
        method_labels{im}, T.E_inf(idx1), T.E_inf(idx2), T.E_inf(idx3));
    fprintf(fid, '& IT & %d & %d & %d \\\\\n', ...
        T.IT(idx1), T.IT(idx2), T.IT(idx3));
    fprintf(fid, '& Halvings & %d & %d & %d \\\\\n', ...
        T.Halvings(idx1), T.Halvings(idx2), T.Halvings(idx3));
end

fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
end
