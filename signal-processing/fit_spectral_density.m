% Tue 20 Jul 22:40:41 CEST 2021
% function [par,S_] = fit_spectral_density(fx,S,par,L,form,model,method)
function [par,S_,resn] = fit_spectral_density(fx,S,par,L,form,model,method,fmax)
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
		fmax = inf;
	end

	S   = cvec(S);
	n   = length(fx);

	df = 1/L;
	fdx = fx>0 & fx < fmax;
	S  = S/sum(S(fdx));
	
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
	S_ = real(S_);
	fdx  = fx>=0;
	fdx0 = fx==0;
	S_ = S_/(df*(0.5*S_(fdx0) + sum(S_(fdx))));

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
		otherwise
			disp(model);
			error('here')
		end
		S_ = real(S_);
		S_ = S_/sum(S_(fdx));
		switch (method)
		case {'lsrel'}
			res = (S(fdx)-S_(fdx))./S_(fdx);
		case {'ls'}
			res = S_(fdx) - S(fdx);
		case {'ll'}
			res    = (log(S_) + S./S_);
			res(1) = 0;
			res    = sum(res(fdx));
		end
	end
end % fit_spectral_density

