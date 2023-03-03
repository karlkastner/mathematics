% Thu 22 Mar 11:11:40 CET 2018
% Karl Kastner, Berlin
%
%% variance of the log-uniform distribution
function s2 = varlogu(a,b)
	% first central moment (mean)
	x1 = meanlogu(a,b);
	% second central moment
	x2 = (b^2-a.^2)/(2*log(b/a));
	% variance
	s2 = x2 - x1.^2;
end

