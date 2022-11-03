% 2021-06-21 15:53:50.278752860 +0200
% radial periodogram
% this is somewhat in between a periodogram for small r and density for large r
%function [Smu,fri,se,count] = periodogram_radial(Shat,L)
% Smu : S   = 1/k int Shat dtheta
% rS  : k*S = int Shat dtheta = k S
function [Smu, Sn, kS, fri, se, count] = periodogram_radial(Shat,L)
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
	frmax = hypot(max(fx),max(fy));
	fri = (0:dfi:ceil(frmax+dfi))';
%	fri  = (0: 

	% this is a linear interpolation
	% each value is proportionally split between the next lower and next larger bin
	% integer part
	rat = fr/dfi;
	bin = floor(rat);
	p   = 1-(rat-bin);
	% 0 is bin 1
	bin = bin+1;

	s = [length(fri),1];
	% half-sum
	kS = 0.5*( accumarray(bin(:),p(:).*Shat(:),s,@sum) ...
	         + accumarray(bin(:)+1,(1-p(:)).*Shat(:), s, @sum) );
	% number of bins
	count = ( accumarray(bin(:),p(:),s,@sum) ...
	        + accumarray(bin(:)+1,(1-p(:)), s, @sum) );
	%n  = accumarray(flat(r),ones(numel(f),1),[],@sum);
	Smu = kS./count;
	Smu(count==0) = 0;
	% normalize
	Sn = Smu./(sum(mid(Smu)).*dfi);

	% TODO interpolate for sd to
	%se = accumarray(flat(r),flat(f),[],@(x) std(x)/sqrt(length(x)));
end % periodogram_radial

