% Tue 20 Jul 23:47:34 CEST 2021
%
% analyze vegetation pattern, either remotely sensed or model generated
%
function analyze(obj)
	stat   = obj.stat;
	S      = obj.S; 
	R      = obj.R;

	b      = cvec(obj.b);

	n      = length(b);
	x      = obj.x();
	fx     = obj.fx();
	fdx    = find(fx>=0);

	switch (obj.opt.normalize)
	case {0,false}
		% nothing to do
	case {1,true}
		b      = b-min(b);
		b      = b/rms(b);
	case {'detrend'}
		x  = (1:n)';
		A  = [ones(n,1),x];
		b_ = A \ b;
		b  = b - A*b_;
		b  = b-min(b);
		b  = b/rms(b);
	end
	b0     = b-mean(b);

	% moments of the raw periodogram
	S.raw0      = periodogram(b0,obj.L);
	f_mean0     = wmean(S.raw0(fdx),fx(fdx));

	fmin = obj.opt.fminscale*f_mean0;
	fmax = obj.opt.fmaxscale*f_mean0;
	% remove low requent lobe through highpass filtering
	if (~isempty(fmin) && obj.fmin < inf)
		p_hp  = 2;
		b_    = highpass1d_fft_cos(b0,fmin,obj.L,p_hp);
		obj.b_hp = b_;
	else
		b_ = b0;
	end

	% periodogram
	S.raw  = periodogram(b_,obj.L);

	% moments
	f_mean   = wmean(S.raw(fdx),fx(fdx));
	f_std	 = wstd(S.raw(fdx),fx(fdx));

	% determine number of segments for bartlett's estimate
	if (~isempty(obj.opt.nb))
		stat.m  = obj.opt.nb;
	else
		mseg    = sqrt(obj.opt.mseg_scale*obj.L*f_mean);
		mseg    = min(mseg,0.5*obj.L*f_mean);
		stat.m  = round(mseg);
	end	
	stat.m_ = round(sqrt(0.5*obj.L*f_mean));
%	Lw_rect = stat.m/obj.L;
%	fcut_rect = rectwin_cutoff_frequency(Lw_rect);
%	Lw_gauss = fcut2Lw_gausswin(fcut_rect);

	% non-parameteric densities
	S.raw0        = periodogram(b0,obj.L);
	S.bartlett0   = periodogram_bartlett(b0,obj.L,stat.m,n,[],false);
	S.bartlett    = periodogram_bartlett(b_,obj.L,stat.m,n,[],false);
	S.gauss       = periodogram_gauss(b_,obj.L,stat.m); %Lw_gauss);
	S.gauss_       = periodogram_gauss(b_,obj.L,stat.m_); %Lw_gauss);
	S.rw          = periodogram_rectwin(b_,obj.L,stat.m); %Lw_gauss);
	S.filt        = ifftshift(meanfilt1(fftshift(S.raw),stat.m));
	S.flat        = spectral_density_flat(obj.L,n);
%obj.L
%stat.m
%clf
%plot([S.filt,S.bartlett,S.gauss,S.rw])
%legend('filt','ba','g','rw')
%pause

	f_C = {'gauss','gauss_','bartlett','filt','raw'};

	for idx=1:length(f_C)
		[stat.Sc.(f_C{idx}),mdx] = max(S.(f_C{idx})(fdx));
		stat.fc.(f_C{idx})       = fx(fdx(mdx));
	end
%	stat.fc.bartlett       = fx(fdx(mdx));
%	[stat.Sc.filt,mdx]     = max(S.filt(fdx));
%	stat.fc.filt           = fx(fdx(mdx));

	% reference density for confidence interval
	if (nargin()>5 && ~isempty(Sref))
		S.ref = Sref;
	else
		S.ref = S.bartlett;
	end
	[stat.Sc.ref,mdx] = max(S.ref);
	stat.fc.ref    = fx(mdx);

	% fit log normal density by method of moments
	[lmu,lsd]         = logn_moment2param(f_mean,f_std);
	par0.lognormal    = [lmu,lsd];

	template = 'gauss_';

	% the small perturbation 1/L is necessary to avoid exceptions when fc = 0
	fc0 = stat.fc.(template) + 1./obj.L;


	% fitted parametric spectral densities
	par0.bandpass = [fc0, ...
			    spectral_density_bandpass_continuous_max2par(fc0,stat.Sc.(template))];
	if (~isfinite(par0.bandpass(2)))
		par0.bandpass(2) = 10;
	end
	[par0.brownian(1),par0.brownian(2)] = spectral_density_brownian_phase_mode2par(fc0,stat.Sc.(template));
	par0.lorentzian = [fc0, spectral_density_lorentzian_max2par(fc0,stat.Sc.(template))];
	if (~isfinite(par0.lorentzian(2)))
		par0.lorentzian(2) = 10;
	end

	f_C     = {'bandpass','brownian','lognormal','lorentzian'};
	model_C = {'bandpass-continuous','brownian-phase','lognormal','lorentzian'};
	for idx=1:length(f_C)
		%[par, S.(f_C{idx})] = fit_spectral_density(fx,S.(obj.opt.template), ...
		%	par0.(f_C{idx}),obj.L,'f',model_C{idx},obj.opt.fitmethod,fmin,fmax);
		fitmethod = 'ls'; 
		[par, S.(f_C{idx})] = fit_spectral_density(fx,S.(template), ...
			par0.(f_C{idx}),obj.L,'f',model_C{idx},fitmethod,fmin,fmax); 
		%[par, S.(f_C{idx})] = fit_spectral_density(fx,S.raw, ...
		%	par,obj.L,'f',model_C{idx},fitmethod,fmin,fmax); 
		[Sc,mdx] = max(S.(f_C{idx}));
		stat.Sc.(f_C{idx}) = Sc;
		% the densities are symmetric at zero hence fit is independent of sign
		stat.fc.(f_C{idx}) = abs(fx(mdx));
		% esp for lorentzian, the paramenter 2 is also independent of sign
		stat.par1.(f_C{idx}) = abs(par(1));
		stat.par2.(f_C{idx}) = abs(par(2));
	end
	%S.lognormal       = lognpdf(fx,lmu,lsd);
	%[stat.fc.lognormal,stat.Sc.lognormal] = logn_mode(lmu,lsd);
	% correction of S for truncation of high frequencies
	%scale             = 1./(sum(S.lognormal(fdx))*obj.df);
	%S.lognormal       = scale*S.lognormal;
	%stat.Sc.lognormal = scale*stat.Sc.lognormal;

	% goodness of fit of parametric densities compared to non-parametric densities
	% TODO makes this part of fit_spectral_density
	fdx_ = fx >= fmin & fx < fmax;
	var_ref   = var(S.ref(fdx_));
	for jdx=1:length(f_C)
		% the mse could be corrected (lower) for uncertainty in bartlett's estiamte
		% we do not do this here, as this in extreme cases leads to
		% an R2 > 1
		stat.mse.(f_C{jdx})     = mean((S.(f_C{jdx})(fdx_)-S.ref(fdx_)).^2);
		stat.r2.(f_C{jdx})      = 1 - stat.mse.(f_C{jdx}) ./ var_ref;
		% check for sufficient goodness of fit
		stat.good.(f_C{jdx})    = ( (stat.r2.(f_C{jdx})>=obj.opt.r2_min) ...
					    && (stat.fc.(f_C{jdx}) >= fmin) ...
					  );
	end
	stat.good.bartlett = (stat.fc.bartlett >= fmin);

	S.c = periodogram_confidence_interval(obj.opt.pp,S.ref,1);

	try
		[stat.p_periodic,ratio,maxShat,mdx_,fdx_,S.mue] = periodogram_test_periodicity(fx,S.raw,stat.m,fmin,fmax);
	catch e
		disp(e);
		stat.p_periodic = NaN;
	end
	stat.periodic = (stat.p_periodic<obj.opt.confidence_level);

	% correlogiram and autocorrelation
	R.raw     = autocorr_fft(cvec(b),[],1);
	R.b       = real(ifft(S.bartlett));
	R.b       = R.b/R.b(1);
	R.bp      = real(ifft(S.bandpass));
	R.bp      = R.bp./R.bp;
	R.ref     = real(ifft(S.ref));
	R.ref     = R.ref/R.ref(1);

	% stationarity test
	% kolmcdf is not a standard matlab package and externally provided
	if (~isempty(which('kolmcdf')))
		[stat.p_stationary,stat.D_stationary,pp,ratio,mdx_] = periodogram_test_stationarity(b,stat.m);
	else
		stat.p_stationary = NaN;
		stat.D_stationary = NaN;
	end

	stat.biomass                  = mean(b);
	stat.wavelength.mean          = 1./f_mean;
	f_C = {'gauss','bartlett','bandpass','brownian','filt','lognormal','lorentzian','raw'};
	for idx=1:length(f_C)
		stat.wavelength.(f_C{idx}) = 1./stat.fc.(f_C{idx});
	end

	% moments of frequency distribution f=k/2pi
	stat.f.mean0                = f_mean0;
	stat.f.mean                 = wmean(S.(obj.opt.template)(fdx),fx(fdx)); 
	stat.f.std                  = wstd(S.(obj.opt.template)(fdx),fx(fdx)); 
	stat.f.skew                 = wskew(S.(obj.opt.template)(fdx),fx(fdx)); 
	stat.f.kurt                 = wkurt(S.(obj.opt.template)(fdx),fx(fdx)); 
	stat.acf_decay_rate         = autocorr_decay_rate(x,R.ref);
	% dummies, computed only for models
	stat.celerity               = [];
	stat.corr                   = [];

	% write back
	obj.S = S;
	obj.R = R;
	obj.stat = stat;
end % Spatial_Pattern::analyze

