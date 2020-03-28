% 2012 Feb 26 17:24 MSK
% Karl Kästner, Berlin

% expects points in rows, col1 = x, col2 = y
% assume nucleus at x0 = 0
function f = f_coulomb(q)
%	f = -2.0./sqrt( epsilon^2 + sum( (q - ones(size(q,1),1)*L).^2, 2) );
	f = -1.0./sqrt( sum(q.^2, 2) + eps);
%	f = -1.0./sqrt( sum(q.^2, 1) );
end % function f_coulomb

