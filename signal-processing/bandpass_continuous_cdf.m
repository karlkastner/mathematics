function cdf = cdf_bandpass_continuous(fx,fc,p)
	scale = spectral_density_bandpass_continuous_scale(fc,p);
	cdf = scale*4^p*fx.*(fx/fc).^(2*p).*hypergeom([2*p,p+1/2],p+3/2,-fx.^2/fc.^2)./(2*p+1);
end

