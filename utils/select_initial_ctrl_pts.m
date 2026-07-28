function P0 = select_initial_ctrl_pts(Q, n)
% SELECT_INITIAL_CTRL_PTS  Choose initial control points from data.
%   P0 = select_initial_ctrl_pts(Q, n)
%
%   Strategy (iii) from Wu & Liu (2024) / Deng & Lin (2014) Eq. (23):
%     P0(1,:) = Q(1,:),  P0(n+1,:) = Q(m+1,:)
%     P0(i,:) = Q(floor((m+1)*i/n) + 1, :)  for i = 1,...,n-1  (1-based)
%
%   Input:
%     Q - (m+1)-by-d data points
%     n - max 0-index of control points, so n+1 control points
%
%   Output:
%     P0 - (n+1)-by-d initial control points

[m1, d] = size(Q);   % m+1
m = m1 - 1;

P0 = zeros(n+1, d);
P0(1,:)   = Q(1,:);
P0(n+1,:) = Q(m+1,:);

for i = 1:n-1
    idx = floor(m1 * i / n) + 1;
    idx = min(idx, m1);
    P0(i+1,:) = Q(idx,:);
end
end
