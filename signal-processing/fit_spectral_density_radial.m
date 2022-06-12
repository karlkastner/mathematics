% Wed  4 May 13:58:08 CEST 2022
% Karl Kästner, Berlin
%
%% fit spectral densities
%
% function [par,S_] = fit_spectral_density(fx,S,par,L,form,model,method)
function [par,Sp,resn] = fit_spectral_density_radial(f,S,par,model,method,fmin,fmax)
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

	fdx = f>fmin & f < fmax;

	% initial values
	if (nargin()<3 || isempty(par))
		par(1) = 0.5*sum(f.*S)./sum(S)
		par(2) = 1
	end

	S = S/sum(flat(S));
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
			Sp = spectral_density_bandpass_2d(f,par(1),par(2));
		otherwise
			disp(model);
			error('Spectral density model not available')
		end
		Sp = real(Sp);
		Sp_ = Sp/sum(flat(Sp));
		% for sanity
		switch (method)
		case {'lsrel'}
			res = f.*(Sp_-S)./Sp_;
			res = res(fdx);
		case {'ls'}
			% the weights are 2 pi fr, but 2 pi is a constant
			res = f(fdx).*(Sp_(fdx) - S(fdx));
		case {'ll'}
			res    = (log(Sp) + S./Sp);
			res(1) = 0
			res    = sum(res(fdx));
		end
	end
end % fit_spectral_density

