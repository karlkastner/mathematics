% Sat 30 Apr 09:27:01 CEST 2022
% Karl Kästner, Berlin
%
%% fit spectral densities
%
% function [par,S_] = fit_spectral_density(fx,S,par,L,form,model,method)
function [par,Sp,resn] = fit_spectral_density_2d(fx,S,par,model,method,fmin,fmax)
	if (nargin<4)
		model = 'bandpass';
	end
	if (nargin()<5)
		method = 'ls';
	end
	if (nargin()<6)
		fmin = 0;
	end
	if (nargin()<7)
		fmax = inf;
	end


	% radial frequency
	fr = hypot(fx,fx');

	fdx = fr>fmin & fr < fmax;

	% initial values
	if (nargin()<3 || isempty(par))
		par(1) = 0.5*sum(flat(fr.*S))./sum(flat(S));
		par(2) = 1;
	end

	switch (method)
	case {'ls'}
		opt = struct();
		%opt.Display = 'off';
		[par,resn,res] = lsqnonlin(@(par) resfun(par), par, [], [], opt); 
	case {'ll'}
		[par] = fminsearch(@(par) resfun(par),par); 
		resn = 0;
	otherwise
		error('here')
	end

	function res = resfun(par)
		switch (model)
		case {'bandpass'}
			Sp = spectral_density_bandpass_2d(fr,par(1),par(2));
		otherwise
			disp(model);
			error('Spectral density model not available')
		end
		Sp = real(Sp);
		S = S/sum(flat(S));
		Sp_ = Sp/sum(flat(Sp));
		% for sanity
		switch (method)
		case {'lsrel'}
			res = (S-Sp_)./Sp_;
			res = res(fdx);
		case {'ls'}
			res = Sp_(fdx) - S(fdx);
		case {'ll'}
			res    = (log(Sp) + S./Sp);
			res(1) = 0;
			res    = sum(res(fdx));
		end
	end
end % fit_spectral_density

