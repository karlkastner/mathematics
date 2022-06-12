% Tue 20 Jul 22:40:41 CEST 2021
% Karl Kästner, Berlin
%
%% fit spectral densities
%
% function [par,S_] = fit_spectral_density(fx,S,par,L,form,model,method)
function [par,S_,resn] = fit_spectral_density(fx,S,par,L,form,model,method,fmin,fmax)
	if (nargin()<3 || isempty(par))
		par = [0.6,1]; 
	end
	if (nargin<5)
		form = 'f';
	end
	if (nargin()<7)
		method = 'll';
	end
	if (nargin()<8)
		fmin = 0;
	end
	if (nargin()<9)
		fmax = inf;
	end

	S   = cvec(S);
	n   = length(fx);

	fdx = fx>fmin & fx < fmax;
	IS = spectral_density_area(fx,S);
	S   = S./IS;
	Sflat = 2*L/n;
	
	switch (method)
	case {'ls','lsrel'}
		opt.Display = 'off';
		[par,resn,res] = lsqnonlin(@(par) resfun(par),par,[],[],opt); 
	case {'ll'}
		[par] = fminsearch(@(par) resfun(par),par); 
		resn = 0;
	otherwise
		error('here')
	end
	%S_   = real(S_);
	%df  = 1/L;
	%fdx  = fx>=0;
	%fdx0 = fx==0;
	%IS = spectral_density_area(fx,S_)
	%S_ = S_/IS;
	%(df*(0.5*S_(fdx0) + sum(S_(fdx))));

	function res = resfun(par)
		switch (model)
		case {'bandpass-continuous'}
			S_ = spectral_density_bandpass_continuous(fx,par(1),par(2));
		case {'bandpass-discrete'}
			dx = L/n;
			S_ = spectral_density_bandpass_discrete(fx,par(1),par(2),dx,form);
		case {'lorentzian'}
			S_ = spectral_density_lorentzian(fx,par(1),par(2));
		case {'brownian-phase'}
			S_ = spectral_density_brownian_phase(fx,par(1),par(2));
		case {'lognormal'}
			S_ = lognpdf(abs(fx),par(1),par(2));
		otherwise
			disp(model);
			error('here')
		end
		S_ = real(S_);
		IS_ = spectral_density_area(fx,S);
		S_ = S_/IS_;
		%sum(S_(fdx));
		switch (method)
		case {'lsrel'}
			res   = (S(fdx)-S_(fdx))./(S_(fdx)+Sflat);
		case {'ls'}
			res = S_(fdx) - S(fdx);
		case {'ll'}
			res = log(S./(S_+Sflat));
			%res    = (log(S_+Sflat) + S./(S_+Sflat));
			res(1) = 0;
			res    = sum(res(fdx));
		end
	end
end % fit_spectral_density

