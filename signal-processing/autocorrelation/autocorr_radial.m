% 2021-06-23 16:22:53.456353906 +0200
% Fri 15 Jul 17:45:22 CEST 2022
function [Rr,ri,count] = autocorr_radial(R,L)
	n = size(R);
	if (nargin()<2)
		L = n;
	end
	x = fftshift(linspace(-L(1)/2,L(1)/2,n(1)));
	y = fftshift(linspace(-L(2)/2,L(2)/2,n(2)))';
	
	r    = hypot(x,y);
	rmax = hypot(L(1)/2,L(2)/2); 
	dri   = hypot(L(1)/n(1),L(2)/n(2));
	ri   = (0:dri:ceil(rmax+dri))';

	% integer part
	rat = r/dri;
	bin = floor(rat);
	p   = 1-(rat-bin);
	% 0 is bin 1
	bin = bin+1;

	s = [length(ri),1];
	% sum
	sum_R = ( accumarray(bin(:),p(:).*R(:),s,@sum) ...
	         + accumarray(bin(:)+1,(1-p(:)).*R(:), s, @sum) );
	% number of bins
	count = ( accumarray(bin(:),p(:),s,@sum) ...
	        + accumarray(bin(:)+1,(1-p(:)), s, @sum) );
	%n  = accumarray(flat(r),ones(numel(f),1),[],@sum);
	Rr = sum_R./count;
	Rr = Rr/Rr(1);

	% TODO interpolate for sd to
	%se = accumarray(flat(r),flat(f),[],@(x) std(x)/sqrt(length(x)));
end % periodogram_radial
