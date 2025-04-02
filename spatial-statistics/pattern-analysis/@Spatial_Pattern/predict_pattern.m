function predict_pattern(obj)
	T = interp1(obj.f.r,obj.T.radial,obj.f.rr,'linear',0);
	obj.b_.lin = ifft2(T.*fft2(obj.source));
	% strip spurious imaginary part introduced by finite machine precision
	obj.b_.lin = real(obj.b_.lin);
	% threshhold
	q = quantile(obj.b_.lin,1-obj.stat.coverage,'all');
	obj.b_.lin_thresh = obj.b_.lin > q;
	% goodness of fit
	obj.stat.fit.b_.lin.r2        = corr(obj.b_.lin(:),obj.b_.square(:)).^2;
	obj.stat.fit.b_.lin_thresh.r2 = sign_to_pearson(corr(obj.b_.lin_thresh(:),obj.b_.thresh(:))).^2;
end

