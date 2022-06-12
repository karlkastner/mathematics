% 2021-06-21 15:53:50.278752860 +0200
% radial periodogram
% this is somewhat in between a periodogram for small r and density for large r
function [Smu,fri,se,count] = periodogram_radial(Shat,L)
	n = size(Shat);
	if (nargin()<2)
		L = n;
	end
	fx = fourier_axis(L(1),n(1));
	fy = fourier_axis(L(2),n(2));
	fmax = max(max(fx),max(fy)); 
	df   = 1./L;
	dfi  = sqrt(df(1).*df(2));
	fr   = hypot(fx,fy');
	fri  = (0:dfi:ceil(fmax))';

	% this is a linear interpolation
	% each value is proportionally split between the next lower and next larger bin
	% integer part
	rat = fr/dfi;
	bin = floor(rat);
	p   = 1-(rat-bin);
	% 0 is bin 1
	bin = bin+1;

	% sum
	s = [length(fri),1];
	sum_ = ( accumarray(bin(:),p(:).*Shat(:),s,@sum) ...
	       + accumarray(bin(:)+1,(1-p(:)).*Shat(:), s, @sum) );
	count = ( accumarray(bin(:),p(:),s,@sum) ...
	       + accumarray(bin(:)+1,(1-p(:)), s, @sum) );
	%n  = accumarray(flat(r),ones(numel(f),1),[],@sum);
	Smu = 0.5*sum_./count;
	Smu(count==0) =0;
	% TODO interpolate for sd to
	%se = accumarray(flat(r),flat(f),[],@(x) std(x)/sqrt(length(x)));
end % periodogram_radial

