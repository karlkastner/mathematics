% 2021-06-21 22:52:10.758702394 +0200 t_spectrum.m
function [mu,se,r,n] = periodogram_annular(f,n,x)
	if (nargin()<3)
		x = (1:length(f));
	end

	x = x-mean(x);
	t = atan(x'./x)+pi;
	t = atan2(x',x)+pi;
 	%r = sqrt(x.^2 + (x.^2)');
 	%nf = 1;
	%r = round(r/nf)*nf+1;
	t = floor(n*t/(2*pi))+1;
	mu = accumarray(flat(t),flat(f),[],@(x) mean(x));
	se = accumarray(flat(t),flat(f),[],@(x) std(x)/sqrt(length(x)));
%	n = accumarray(flat(t),ones(numel(f),1),[],@sum);
end

