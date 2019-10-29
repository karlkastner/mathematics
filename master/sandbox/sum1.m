% Mon Oct 31 21:01:35 MSK 2011
% Karl Käsnter, Berlin

% sum formula for single sum (gauss)
% sum i
function s = sum1(n)
	p = [1 1 0]/2;
	s = polyval(p,n);
end


