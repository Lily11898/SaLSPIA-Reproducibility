%% ========================================================================
%  run_table7_missing_surface.m
%  Rank-deficient (missing-data) SURFACE fitting on the Franke surface.
%
%  Pipeline:
%    1. sample the clean Franke surface (no noise);
%    2. punch a few local elliptical holes in the (u,v) domain, so some
%       control points lose data support -> the observed collocation matrix
%       is rank-deficient (CONFIRMED numerically via rank(A_obs) < n_ctrl);
%    3. fit from the OBSERVED samples only with four methods
%       (LSPIA / LSPIA-Lin2018 / ALSPIA / SaLSPIA);
%    4. report held-out hole errors (Hole RMS, Hole Max) + full-grid RMSE,
%       plus IT, CPU (averaged over repeats), E_inf;
%    5. draw a convergence figure and a two-panel result figure
%       (data+holes | fitted surface) with numbered zoom insets.
%
%  Fairness notes:
%    - Initial control net uses OBSERVED data only: x,y are fixed grid
%      coordinates; the missing z is filled by interpolation/extrapolation
%      from observed samples, so no held-out hole-interior truth leaks in.
%    - CPU is a single-machine timing; IT is the more reliable speed metric.
%    - A_obs is dense on purpose (solvers call eig(A'*A), which rejects sparse).
%
%  Run setup_paths.m once before this script.
%% ========================================================================
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
repo_dir   = fileparts(script_dir);
addpath(fullfile(repo_dir, 'algorithms'));
addpath(fullfile(repo_dir, 'utils'));

%% ===== Config ==========================================================
save_results = true;       % write CSV + MAT + LaTeX
save_figure  = false;      % write figures (PNG/FIG)
out_dir      = fullfile(repo_dir, 'results', 'missing_data_surface');

p_deg   = 3;   tol_Ek = 1e-6;   maxiter = 10000;
m1u = 101;  m1v = 101;   n1u = 41;  n1v = 41;      % data grid / control mesh
nrep = 5;                                          % CPU repeats (converged only)

params.c = 1e-4;  params.M = 10;  params.delta = 1e-8;   % SaLSPIA defaults
params.eps_saf = 1e-30;                                  % Algorithm 1 safeguard

% Local holes in the (u,v) domain: [cu cv su sv]
% Placed on the gentler slopes of Franke, away from the sharp peaks/pit,
% so the missing regions are reconstructable rather than pure extrapolation.
hole_specs = [ 0.45 0.25 0.11 0.12;
               0.72 0.62 0.11 0.11;
               0.25 0.55 0.12 0.11;
               0.60 0.85 0.11 0.11 ];
SX = 1;  SY = 1;                                   % u->x, v->y scales (Franke: x=u, y=v)

%% ===== Build surface, holes, observed system ===========================
fprintf('\n============================================================\n');
fprintf('Rank-deficient surface fitting: Franke with local holes\n');
fprintf('============================================================\n');

[Qx, Qy, Qz, u_par, v_par] = generate_franke_surface(m1u, m1v);
knots_u = make_clamped_knots(n1u, p_deg);
knots_v = make_clamped_knots(n1v, p_deg);
Au = build_collocation(knots_u, p_deg, u_par);
Av = build_collocation(knots_v, p_deg, v_par);

[U, V] = ndgrid(u_par, v_par);
hole_mask = false(size(U));
for h = 1:size(hole_specs,1)
    c = hole_specs(h,:);
    in = ((U-c(1))/c(3)).^2 + ((V-c(2))/c(4)).^2 <= 1;
    hole_mask = hole_mask | in;
    fprintf('  hole %d: center (u,v)=(%.2f,%.2f), removed %d pts\n', ...
        h, c(1), c(2), nnz(in));
end
obs_mask = ~hole_mask;

A_obs = build_observed_collocation(Au, Av, obs_mask);
Q_obs = [Qx(obs_mask), Qy(obs_mask), Qz(obs_mask)];

% Hole-blind initial control net (observed-only fill of the missing z).
Qz_fill = fill_from_observed(Qz, obs_mask, Qx, Qy);
P0 = [select_initial_surf(Qx,      n1u,n1v), ...
      select_initial_surf(Qy,      n1u,n1v), ...
      select_initial_surf(Qz_fill, n1u,n1v)];

n_ctrl = n1u*n1v;  rankA = rank(A_obs);
n_dead = nnz(sum(A_obs,1) <= 1e-12*max(1,max(sum(A_obs,1))));
fprintf('  Grid %dx%d, controls %dx%d (%d)\n', m1u,m1v,n1u,n1v,n_ctrl);
fprintf('  observed %d, missing %d\n', nnz(obs_mask), nnz(hole_mask));
fprintf('  rank(A_obs) = %d / %d  (inactive columns: %d)  ->  %s\n', ...
    rankA, n_ctrl, n_dead, ternary(rankA<n_ctrl,'RANK-DEFICIENT','full rank'));

%% ===== Run the four methods (CPU averaged over nrep for converged) ======
solvers = {
    'LSPIA',         @() lspia(A_obs, Q_obs, P0, tol_Ek, maxiter);
    'LSPIA-Lin2018', @() lspia_lin2018(A_obs, Q_obs, P0, tol_Ek, maxiter);
    'ALSPIA',        @() alspia(A_obs, Q_obs, P0, tol_Ek, maxiter);
    'SaLSPIA',       @() salspia(A_obs, Q_obs, P0, tol_Ek, maxiter, params);
};
nM = size(solvers,1);
res   = struct('name',{},'IT',{},'CPU_mean',{},'CPU_std',{}, ...
    'Einf',{},'HoleRMS',{},'HoleMax',{},'RMSE',{},'P',{});
infos = cell(1,nM);
cpu_samples_all = cell(1,nM);

for i = 1:nM
    fprintf('  running %-14s ...', solvers{i,1});
    [P, info] = solvers{i,2}();                 % IT / E_inf / P (deterministic)
    cpu_samples = info.cpu_time;
    if info.iter_count < maxiter && nrep > 1     % average CPU only if converged
        cpu_samples = zeros(nrep,1);
        cpu_samples(1) = info.cpu_time;
        for ir = 2:nrep
            [~, info_ir] = solvers{i,2}();
            cpu_samples(ir) = info_ir.cpu_time;
        end
        cpu_mean = mean(cpu_samples);
        cpu_std = std(cpu_samples, 0);
    else
        cpu_mean = NaN;
        cpu_std = NaN;
    end
    info.cpu_time_mean = cpu_mean;
    info.cpu_time = cpu_mean;
    info.cpu_time_std = cpu_std;
    info.cpu_time_samples = cpu_samples;
    info.cpu_n_repeats = numel(cpu_samples);
    [hrms,hmax,rmse] = surface_metrics(P, Qx,Qy,Qz, Au,Av, n1u,n1v, hole_mask);
    res(i)  = struct('name',solvers{i,1},'IT',info.iter_count, ...
        'CPU_mean',cpu_mean,'CPU_std',cpu_std,'Einf',info.Ek_history(end), ...
        'HoleRMS',hrms,'HoleMax',hmax,'RMSE',rmse,'P',P);
    infos{i} = info;
    cpu_samples_all{i} = cpu_samples;
    if info.iter_count >= maxiter
        fprintf(' FAILED (IT>=%d)\n', maxiter);
    else
        fprintf(' IT=%d  CPU=%.6f +/- %.6f  E=%.2e\n', ...
            info.iter_count, cpu_mean, cpu_std, info.Ek_history(end));
    end
end

%% ===== Summary table ===================================================
fprintf('\n%s\n', repmat('=',1,116));
fprintf('  SUMMARY: Rank-Deficient Surface Fitting (Franke)   rank = %d / %d\n', rankA, n_ctrl);
fprintf('%s\n', repmat('=',1,116));
fprintf('%-16s | %8s | %12s | %12s | %10s | %11s | %11s | %11s\n', ...
    'Method','IT','CPU mean(s)','CPU SD(s)','E_inf','Hole RMS','Hole Max','RMSE');
fprintf('%s\n', repmat('-',1,116));
for i = 1:nM
    r = res(i);
    if r.IT >= maxiter
        fprintf('%-16s | %8s | %12s | %12s | %10s | %11s | %11s | %11s\n', ...
            r.name,'#','#','#','#','-','-','-');
    else
        fprintf('%-16s | %8d | %12.6f | %12.6f | %10.2e | %11.4e | %11.4e | %11.4e\n', ...
            r.name, r.IT, r.CPU_mean, r.CPU_std, r.Einf, ...
            r.HoleRMS, r.HoleMax, r.RMSE);
    end
end
fprintf('%s\n', repmat('=',1,116));

if save_results
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    T = struct2table(rmfield(res,'P'));
    writetable(T, fullfile(out_dir,'missing_surface_summary.csv'));
    write_franke_latex(T, fullfile(out_dir,'table7_franke_cpu.tex'), maxiter);
    save(fullfile(out_dir,'missing_surface_results.mat'), 'res','T', ...
        'Qx','Qy','Qz','obs_mask','hole_mask','Au','Av','u_par','v_par', ...
        'knots_u','knots_v','hole_specs','rankA','params','cpu_samples_all');
    fprintf('saved CSV + MAT + LaTeX to: %s\n', out_dir);
end

%% ===== Figures =========================================================
plot_convergence(res, infos, tol_Ek, rankA, n_ctrl, out_dir, save_figure);
plot_result(Qx,Qy,Qz, obs_mask,hole_mask, res(4).P, Au,Av, n1u,n1v, ...
    hole_specs, SX,SY, out_dir, save_figure);


%% ======================= helper functions ==============================
function [Qx,Qy,Qz,u_par,v_par] = generate_franke_surface(m1u,m1v)
% Franke's function on [0,1]x[0,1] (a classic surface-fitting test height field).
u = linspace(0,1,m1u);  v = linspace(0,1,m1v);
[T1,T2] = meshgrid(u,v);  T1 = T1';  T2 = T2';
Qz = 0.75*exp(-((9*T1-2).^2 + (9*T2-2).^2)/4) ...
   + 0.75*exp(-((9*T1+1).^2)/49 - (9*T2+1)/10) ...
   + 0.50*exp(-((9*T1-7).^2 + (9*T2-3).^2)/4) ...
   - 0.20*exp(-(9*T1-4).^2 - (9*T2-7).^2);
Qx = T1;  Qy = T2;
u_par = (u(:)-u(1))/(u(end)-u(1));  v_par = (v(:)-v(1))/(v(end)-v(1));
end

function knots = make_clamped_knots(n1,p)
n_int = (n1-1)-p;  knots = zeros(1,n1+p+1);  knots(end-p:end) = 1;
if n_int > 0, in = linspace(0,1,n_int+2);  knots(p+2:end-p-1) = in(2:end-1); end
end

function A = build_observed_collocation(Au,Av,obs_mask)
[iu,iv] = find(obs_mask);  A = zeros(numel(iu), size(Au,2)*size(Av,2));
for k = 1:numel(iu), A(k,:) = kron(Av(iv(k),:), Au(iu(k),:)); end
end

function Vf = fill_from_observed(Vg, obs_mask, Qx, Qy)
F = scatteredInterpolant(Qx(obs_mask), Qy(obs_mask), Vg(obs_mask), 'linear','nearest');
Vf = reshape(F(Qx(:),Qy(:)), size(Vg));
end

function c = select_initial_surf(Qc,n1u,n1v)
iu = round(linspace(1,size(Qc,1),n1u));  iv = round(linspace(1,size(Qc,2),n1v));
c  = reshape(Qc(iu,iv), [], 1);
end

function [Sx,Sy,Sz] = eval_surface(P,Au,Av,n1u,n1v)
Sx = Au*reshape(P(:,1),n1u,n1v)*Av';
Sy = Au*reshape(P(:,2),n1u,n1v)*Av';
Sz = Au*reshape(P(:,3),n1u,n1v)*Av';
end

function [hrms,hmax,rmse] = surface_metrics(P,Qx,Qy,Qz,Au,Av,n1u,n1v,hole_mask)
[Sx,Sy,Sz] = eval_surface(P,Au,Av,n1u,n1v);
d2 = (Sx-Qx).^2 + (Sy-Qy).^2 + (Sz-Qz).^2;
rmse = sqrt(mean(d2(:)));  dh = sqrt(d2(hole_mask));
hrms = sqrt(mean(dh.^2));  hmax = max(dh);
end

function out = ternary(cond,a,b), if cond, out=a; else, out=b; end, end

function safe_export(fig,fname)
try
    exportgraphics(fig,fname,'Resolution',300);
catch
    saveas(fig,fname);
end
end

%% ---- convergence figure ----------------------------------------------
function plot_convergence(res, infos, tol_Ek, rankA, n_ctrl, out_dir, save_figure)
sty = {'-o','-d','-s','-p'};
col = [0.20 0.42 0.70; 0.55 0.38 0.66; 0.28 0.66 0.48; 0.87 0.36 0.11];
fig = figure('Color','w','Position',[100 100 720 520]); hold on;
for i = 1:numel(infos)
    y = infos{i}.Ek_history;  y(~isfinite(y)|y<=0) = realmin;
    it = infos{i}.iter_count;
    mk = 1:max(1,floor(it/14)):it+1;
    semilogy(0:it, y, sty{i}, 'Color',col(i,:), 'MarkerFaceColor',col(i,:), ...
        'MarkerSize',5+2*(i==4), 'LineWidth',1.1+0.7*(i==4), ...
        'MarkerIndices',mk, 'DisplayName',res(i).name);
end
set(gca,'YScale','log');
yline(tol_Ek,'--','Color',[0.55 0.55 0.55],'HandleVisibility','off');
hold off; grid on; set(gca,'FontName','Helvetica','FontSize',11,'GridAlpha',0.12);
xlabel('Iteration $k$','Interpreter','latex','FontSize',12);
ylabel('$E_k$','Interpreter','latex','FontSize',12);
title(sprintf('Rank-deficient surface (rank = %d / %d)', rankA, n_ctrl), 'FontWeight','normal');
legend('Location','northeast','Box','off','FontSize',10);
if save_figure
    base = fullfile(out_dir,'missing_surface_convergence');
    safe_export(fig,[base '.png']); savefig(fig,[base '.fig']);
end
end

%% ---- result figure: (a) data+holes | (b) fitted surface, + insets -----
function plot_result(Qx,Qy,Qz, obs_mask,hole_mask, P, Au,Av, n1u,n1v, ...
    hole_specs, SX,SY, out_dir, save_figure)

nH = size(hole_specs,1);  VIEW = [42 28];  th = linspace(0,2*pi,200);
[Sx,Sy,Sz] = eval_surface(P,Au,Av,n1u,n1v);
Fz = griddedInterpolant({Qx(:,1), Qy(1,:).'}, Sz, 'linear','nearest');

col_obs=[0.42 0.58 0.78]; col_miss=[0.86 0.18 0.22]; col_lbl=[0.94 0.52 0.14];

% per-hole geometry: center, zoom window, ring (on fitted surface), z-window
H = struct('cx',{},'cy',{},'cz',{},'rx',{},'ry',{},'zl',{},'rgx',{},'rgy',{},'rgz',{});
for h = 1:nH
    c  = hole_specs(h,:);
    cx = SX*c(1);  cy = SY*c(2);
    rx = 1.7*SX*c(3); ry = 1.7*SY*c(4);
    rgx = cx + SX*c(3)*cos(th);  rgy = cy + SY*c(4)*sin(th);
    win = abs(Qx-cx)<=rx & abs(Qy-cy)<=ry;
    zl  = [min(Sz(win)) max(Sz(win))];  if diff(zl)<1e-6, zl=zl+[-0.5 0.5]; end
    pad = 0.12*diff(zl);
    H(h) = struct('cx',cx,'cy',cy,'cz',Fz(cx,cy),'rx',rx,'ry',ry, ...
        'zl',[zl(1)-pad zl(2)+pad],'rgx',rgx,'rgy',rgy,'rgz',Fz(rgx,rgy)+0.004);
end

fig = figure('Color','w','Position',[60 80 1460 600]);
mw=0.30; iw=0.095; gap=(1-2*mw-2*iw)/5; mb=0.14; mh=0.78;
xa=gap; xia=xa+mw+gap; xb=xia+iw+gap; xib=xb+mw+gap;
ig=0.014; ih=min(0.165,(mh-(nH-1)*ig)/nH);

% (a) data with holes
axa = axes('Position',[xa mb mw mh]);
ho = scatter3(Qx(obs_mask),Qy(obs_mask),Qz(obs_mask),4,col_obs,'filled', ...
    'MarkerFaceAlpha',0.55,'MarkerEdgeColor','none'); hold on;
hm = scatter3(Qx(hole_mask),Qy(hole_mask),Qz(hole_mask),9,col_miss,'o','LineWidth',0.6);
for h = 1:nH
    text(H(h).cx,H(h).cy,H(h).cz+0.12,sprintf('%d',h),'Color','w','FontWeight','bold', ...
        'FontSize',10,'HorizontalAlignment','center','BackgroundColor',col_lbl,'Margin',2);
end
hold off; style3(axa,VIEW,true);
legend(axa,[ho hm],{'Observed data','Missing region'},'Location','northeast','Box','off','FontSize',8);
xl=axa.XLim; yl=axa.YLim; zl=axa.ZLim;
for h = 1:nH
    axi = axes('Position',[xia mb+mh-h*(ih+ig) iw ih]);
    scatter3(Qx(obs_mask),Qy(obs_mask),Qz(obs_mask),3,col_obs,'filled', ...
        'MarkerFaceAlpha',0.55,'MarkerEdgeColor','none'); hold on;
    scatter3(Qx(hole_mask),Qy(hole_mask),Qz(hole_mask),7,col_miss,'o','LineWidth',0.6); hold off;
    xlim([H(h).cx-H(h).rx H(h).cx+H(h).rx]); ylim([H(h).cy-H(h).ry H(h).cy+H(h).ry]); zlim(H(h).zl);
    style3(axi,VIEW,false); badge(h,col_lbl);
end

% (b) fitted surface + hole rings (white halo so red reads on jet)
axb = axes('Position',[xb mb mw mh]);
hf = surf(Sx,Sy,Sz,'EdgeColor','none','FaceColor','interp','FaceLighting','gouraud', ...
    'AmbientStrength',0.6,'DiffuseStrength',0.6,'SpecularStrength',0.08); hold on;
colormap(axb,jet); camlight(-32,62);
hr = gobjects(1);
for h = 1:nH
    plot3(H(h).rgx,H(h).rgy,H(h).rgz,'-','Color','w','LineWidth',3);
    hln = plot3(H(h).rgx,H(h).rgy,H(h).rgz,'-','Color',col_miss,'LineWidth',1.5);
    if h==1, hr=hln; end
end
hold off; style3(axb,VIEW,false); xlim(xl); ylim(yl); zlim(zl);
legend(axb,[hr hf],{'Missing region','Fitted surface'},'Location','northeast','Box','off','FontSize',8);
for h = 1:nH
    axi = axes('Position',[xib mb+mh-h*(ih+ig) iw ih]);
    surf(Sx,Sy,Sz,'EdgeColor','none','FaceColor','interp','FaceLighting','gouraud', ...
        'AmbientStrength',0.6,'DiffuseStrength',0.6,'SpecularStrength',0.08); hold on;
    colormap(axi,jet); camlight(-32,62);
    plot3(H(h).rgx,H(h).rgy,H(h).rgz,'-','Color','w','LineWidth',2.2);
    plot3(H(h).rgx,H(h).rgy,H(h).rgz,'-','Color',col_miss,'LineWidth',1.2); hold off;
    xlim([H(h).cx-H(h).rx H(h).cx+H(h).rx]); ylim([H(h).cy-H(h).ry H(h).cy+H(h).ry]); zlim(H(h).zl);
    style3(axi,VIEW,false); badge(h,col_lbl);
end

if save_figure
    base = fullfile(out_dir,'missing_surface_fit');
    safe_export(fig,[base '.png']); savefig(fig,[base '.fig']);
end
end

function style3(ax,VIEW,mainpanel)
view(ax,VIEW); grid(ax,'on'); box(ax,'on');
if mainpanel, axis(ax,'tight'); end
set(ax,'FontName','Helvetica','FontSize',7+mainpanel,'GridAlpha',0.12, ...
    'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'ZColor',[.3 .3 .3], ...
    'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
ax.Clipping = ternary(mainpanel,'off','on');
end

function caption(fig,x0,y0,w,str)
annotation(fig,'textbox',[x0 y0-0.085 w 0.05],'String',str,'EdgeColor','none', ...
    'HorizontalAlignment','center','FontName','Helvetica','FontSize',10,'Color',[.15 .15 .15]);
end

function badge(idx,col)
text(0.08,0.92,sprintf('%d',idx),'Units','normalized','Color','w','FontWeight','bold', ...
    'FontSize',8,'BackgroundColor',col,'Margin',2,'VerticalAlignment','top');
end

function write_franke_latex(T, fname, maxiter)
fid = fopen(fname, 'w');
if fid == -1
    warning('Cannot write Franke CPU LaTeX table.');
    return;
end

fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Missing-data results for the Franke surface in Example~4.8. CPU times are the mean $\\pm$ sample standard deviation over five runs.}\n');
fprintf(fid, '\\label{tab:missing_franke_cpu_sd}\n');
fprintf(fid, '\\begin{tabular}{lrrrrr}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Method & IT & CPU (s), mean $\\pm$ SD & $E_\\infty$ & Hole RMS & Hole Max \\\\\n');
fprintf(fid, '\\hline\n');
for i = 1:height(T)
    if T.IT(i) >= maxiter
        fprintf(fid, '%s & $>10^4$ & \\# & \\# & -- & -- \\\\\n', T.name{i});
    else
        fprintf(fid, '%s & %d & $%.6f \\pm %.6f$ & %.2e & %.3e & %.3e \\\\\n', ...
            T.name{i}, T.IT(i), T.CPU_mean(i), T.CPU_std(i), ...
            T.Einf(i), T.HoleRMS(i), T.HoleMax(i));
    end
end
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
end
