% Mon  6 Sep 10:10:13 CEST 2021
% Karl Kästner, Berlin
%
%% lorentzian spectral density
%
% function S = spreactral_density_lorentzian(fx,varargin)
function [S,IS] = spreactral_density_lorentzian(fx,fc,p,normalize)
	if (nargin()<4)
		normalize = true;
	end
%	if (2 == length(varargin))
%		f0 = varargin{1};
%		p  = varargin{2};
%	else
%		f0 = varargin{1}(1);
%		p  = varargin{1}(2);
%	end

	S = 1./(1 + (p.*(abs(fx./fc)-1)).^2);
	% normalize
	switch (normalize)
	case {-1}
		% do not normalize
	case {0,false}
		% normalize analytically
		IS = spectral_density_lorentzian_scale(fc,p);
		S  = S./IS;
	otherwise % numerically
		IS = spectral_density_area(fx,S);
		S  = S./IS;
	end
end
