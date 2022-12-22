% Thu  6 Oct 14:38:53 CEST 2022
% Wed 16 Nov 10:29:36 CET 2022
%
% estimate the maximum of the spectral density by fitting a quadratic polynomial
% this reduces the bias
%
% TODO if S is a radial density, then values should be weighted proportionally to r
%
function [fc,Sc,reg,rmse,Sp] = spectral_density_maximum_bias_corrected(fx,Shat,opt)
	fx   = cvec(fx);
	Shat = cvec(Shat);

	if (nargin()<3)
		opt = struct();
	end
	[ifx,intS] = spectral_density_cumulative_integral(fx,Shat);
	% determine smoothing window width length as range containing 50% of the spectral energy
	if (~isfield('fmin',opt))
		opt.fmin = interp1(ifx,intS,0.25,'linear');
	end
	if (~isfield('fmax',opt))
		opt.fmax = interp1(ifx,intS,0.75,'linear');
	end
	dx = fx(2)-fx(1);
	nf = (opt.fmax - opt.fmin)/dx;

	if (isfiled(opt,'logx') && opt.logx)
		x = log(fx);
	else
		x = fx;
	end
	if (isfield(opt,'logy') && opt.logy)
		y = log(Shat);
	else
		y = Shat;
	end

	% match by variance : gaussian with same standard deviation as rectangular window
	% as sum of squared weights is proportional to dof when filtering uncorrelated data
	% nota bene: this is only approximate, as the the weights are discrete
	%            but the discrepancy is negligible
	sigma = (fx(2)-fx(1))*nf/sqrt(12)
	% complete symmetric half plane, if density is only given for the positive half 
	[fx,Shat] = fourier_complete_negative_half_plane(fx,Shat); 
	n    = length(Shat);
	% initial estimate of the maximum
	% smoothing window
	w    = normpdf(fx,0,sigma);
	% normalize weights to 1
	w    = w/sum(w);
	% density
	Sbar = cvec(ifftshift(conv(fftshift(Shat),fftshift(w),'same')));
	% maximum value
	[Sc, mdx] = max(Sbar);
	% frequency of the maximum
	fc_   = fx(mdx);
	fc_   = abs(fc_);
	fc.biased = fc;
	Sc.biased = Sc;
	reg_biased = Sc.biased.*fc.biased;

	% weights for restriction to the neighbourhood of the maximum
	% (equivivalent to smoothing)
	w = normpdf(fx,fc,sigma);
	% apply radial weights to account for dof dependence on radius
	if (isfield(opt.isradial) && opt.isradial)
		w = abs(fx).*w;
	end
	% normalize sum of weights to 1
	w = w/sum(w);
	W    = diag(sparse(w));

	% re-estimate of the maximum of the density by fitting a quadratic polynomial
	% note: the estimate can be improved by iteration
	% vandermonde matrix of regression
	A    = vander_1d(x,2);

	% coefficients of the quadratic polynomial
	cc   = (A'*W*A) \ (A'*W*y);

	% predicted density (only valid near fc)
	yp = A*cc;

	% residual
	res  = y - yp;

	% mean square error of S
	% normalization not necessary, as sum w = 1
	s2.y = res'*(w.*res);

	% covariance matrix
	C   = s2.y*inv(A'*W*A);

	a = cc(3);
	b = cc(2);
	c = cc(1);

	% frequency fc of the maximum
	xc = -0.5*cc(2)/cc(3);

	% density Sc at the maximum
	% vandermonde matrix at the maximum
	%A    = vander_1d(fc,2);
	%yc   = A*cc;
	yc = c - b^2/(4*a);

	% elements of the covariance matrix
	sa2 = C(3,3);
	sb2 = C(2,2);
	sc2 = C(1,1);
	sab = C(3,2);
	sac = C(3,1);
	sbc = C(2,1);

	s2.lreg = ( 16*a^4*sc2 ...
	-  16*a^3*(b-1)*sbc ...
	+  8*a^2*b*(b-2)*sac ...
	+ sb2*4*a^2*(b^2 - 2*b + 1) ...
	+ sab*4*a*b*(-b^2 + 3*b - 2) ...
	+ sa2*b^2*(b^2 - 4*b + 4) ...
	)/(16*a^4);

	s2.xc = (...
	sb2/(4*a^2) ...
	+ (sa2*b^2)/(4*a^4) ...
	- (sab*b)/(2*a^3) ...
	);

	s2.yc = (sc2 ...
	 + (sa2*b^4)/(16*a^4) ...
	 + sb2*(b^2)/(4*a^2) ...
	 - sbc*(b)/a ...
	 - sab*(b^3)/(4*a^3) ...
	 + sac*(b^2)/(2*a^2) );

	if (opt.logx)
		fc.corrected = exp(xc);
	        rmse.fc      = sqrt(s2.xc).*fc.corrected; 
	else
		fc.corrected = xc;
	        rmse.fc      = sqrt(s2.xc); 
	end

	if (opt.logy)
		Sc.corrected = exp(yc);
		rmse.Sc  = sqrt(s2.yc).*Sc.corrected; 
		Sp       = exp(yp);	
	else
		Sc.corrected = yc;
		rmse.Sc  = sqrt(s2.yc); 
		Sp = yp;
	end

	% regularity
	% note that reg is still biased by error propagation of product
	% TODO this was incorrectly -, check error estimate
	reg.corrected = Sc.corrected*fc.corrected;
	% TODO this is derived for when x and y have been logged
	rmse.reg.corrected = sqrt(s2.lreg)*reg.corrected;
end

