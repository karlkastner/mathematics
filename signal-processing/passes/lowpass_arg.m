% Tue 27 Jul 20:02:12 CEST 2021
% Karl Kästner, Berlin
% maximum frequency for the symmetric lowpass (rho = rho_lp = rho_hp)
% function [r,rho] = lowpass_arg(arg,dx,form,varargin)	
function [r,rho] = lowpass_arg(arg,dx,p,form)
	if (nargin < 4)
		form = 'rho';
	end
	switch (form)
	case {'omega'}
		omega0 = arg;
		fc     = omega/(2*pi);
		rho    = lowpass_fc_to_rho(fc,p,dx);
		r      = rho./(1-2*rho+rho.*rho);
	case {'f'}
		fc     = arg;
		rho    = lowpass_fc_to_rho(fc,p,dx);
		r      = rho./(1-2*rho+rho.*rho);
	case {'rho'}
		rho    = arg;
		r      = rho./(1-2*rho+rho.*rho);
	case {'r'}
		r      = arg;
	otherwise
		error('');
	end
end

