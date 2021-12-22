% Tue 27 Jul 20:02:12 CEST 2021
% function [r,rho] = filter_arg(arg,dx,form,varargin)	
function [r,rho] = filter_arg(arg,dx,form)	
	if (nargin < 3)
		form = 'rho';
	end
	switch (form)
	case {'omega'}
		omega0 = arg;
		f0     = omega/(2*pi);
		rho    = filter_f0_to_rho(f0,dx);
		r      = rho./(1-2*rho+rho.*rho);
	case {'f'}
		f0     = arg;
		rho    = filter_f0_to_rho(f0,dx);
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

