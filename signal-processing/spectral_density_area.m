% 2021-10-21 09:29:25.938016430 +0200
function I = spectral_density_area(fx,S)
	if (isvector(S))
		S = cvec(S);
	end
	df   = fx(2)-fx(1);
	fdx  = fx>0;
	fdx0 = fx==0;
	if (sum(fdx0)>0)
	I    = (0.5*S(fdx0,:)+sum(S(fdx,:)))*df;
	else
	I    = (sum(S(fdx,:)))*df;
	end
end

