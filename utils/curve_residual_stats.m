function stats = curve_residual_stats(A, Q, P)
%CURVE_RESIDUAL_STATS  Geometric residual statistics for curve fitting.
%   stats = CURVE_RESIDUAL_STATS(A,Q,P) computes pointwise residual norms
%   between the fitted values A*P and the data Q.

R = A * P - Q;
pointwise = sqrt(sum(R.^2, 2));

stats.fro_error = norm(R, 'fro');
stats.sse = stats.fro_error^2;
stats.rms_residual = sqrt(mean(pointwise.^2));
stats.max_pointwise_error = max(pointwise);
stats.mean_pointwise_error = mean(pointwise);
stats.pointwise = pointwise;
end
