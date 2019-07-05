% 2015-08-25 12:45:32.053852527 +0200

% stanard error of x with respect to mean

function serr = serr1(x)
	n    = length(x);
	mu   = (1/n)*sum(x);
	d    = x-mu;
	% faster than sum(d.*d)
	s2   = (1/n)*(d'*d);
	serr = sqrt(1/(n-1)*s2);
end

