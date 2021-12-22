% 2021-06-21 15:53:50.278752860 +0200 r_spectrum.m
function [mu,se,r,x,y] = r_spectrum(f,x)
	if (nargin()<2)
		x = (1:size(f,1));
		y = (1:size(f,2));
	end

	x = x-mean(x);
	y = y-mean(y);
 	r = sqrt(x.^2 + (y.^2)')';
 	nf = 1;
	r = round(r/nf)*nf;
	if (mod(size(r,1),2)==1)
		r = r+1;
	end
	mu = accumarray(flat(r),flat(f),[],@(x) mean(x));
	se = accumarray(flat(r),flat(f),[],@(x) std(x)/sqrt(length(x)));
	n  = accumarray(flat(r),ones(numel(f),1),[],@sum);
end


