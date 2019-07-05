% Tue Jul 22 18:49:58 WIB 2014
% Karl Kastner, Berlin
%
%% trimmed mean
function m = qmean(v,q1,q2)
	q = quantile(v, [q1, q2]);
	fdx = find(q(1) < v & v < q(2));
	m = mean(v(fdx));
end

