function stats = surface_residual_stats(Au, Av, Q, P)
%SURFACE_RESIDUAL_STATS  Geometric residual statistics for surface fitting.
%   stats = SURFACE_RESIDUAL_STATS(Au,Av,Q,P) computes pointwise residual
%   norms between tensor-product fitted values Au*P{ell}*Av' and Q{ell}.

d = numel(Q);
err2 = zeros(size(Q{1}));
fro2 = 0;

for ell = 1:d
    R = Au * P{ell} * Av' - Q{ell};
    err2 = err2 + R.^2;
    fro2 = fro2 + norm(R, 'fro')^2;
end

pointwise = sqrt(err2);

stats.fro_error = sqrt(fro2);
stats.sse = fro2;
stats.rms_residual = sqrt(mean(pointwise(:).^2));
stats.max_pointwise_error = max(pointwise(:));
stats.mean_pointwise_error = mean(pointwise(:));
stats.pointwise = pointwise;
end
