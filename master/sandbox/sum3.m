% Mon Oct 31 21:01:14 MSK 2011
% Karl Kästner, Berlin

% sum formula for triple sum
% sum sum sum i
function s = sum3(n)
	p = [1 6 11 6 0] / 24;
	s = polyval(p,n);
end % sum3

