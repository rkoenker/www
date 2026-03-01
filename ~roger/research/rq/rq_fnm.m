function b = rq_fnm(X, y, p)
% Construct the dual problem of quantile regression
% Solve it with lp_fnm
% 
% Function rq_fnm of Daniel Morillo & Roger Koenker
% Found at: http://www.econ.uiuc.edu/~roger/rqn/rq.ox
% Translated from Ox to Matlab by Paul Eilers 1999
%
[m n] = size(X);
u = ones(m, 1);
a = (1 - p) .* u;
b = -lp_fnm(X', -y', X' * a, u, a)';

