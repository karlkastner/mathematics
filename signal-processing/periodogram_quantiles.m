% Wed 28 Jul 17:21:42 CEST 2021
function q = periodogram_quantiles(fx,S,p,smooth)
	if (nargin()<3||isempty(p))
		p = [0.16,0.84];
	end
	if (nargin()>3 && smooth)
		[Smax,mdx] = max(S);
		S = trifilt1(S,2*mdx);
	end

	fdx = fx>0;
	fx = fx(fdx);
	S  = S(fdx);
	Si = cumsum(S);
	Si = Si/Si(end);
	% TODO, smooth
	for idx=2:length(Si); Si(idx)=max(Si(idx),Si(idx-1)*(1+sqrt(eps))); end
	q = interp1(Si,fx,p,'linear');
end


