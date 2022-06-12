% Fri 22 Apr 13:28:53 CEST 2022
% fr : radial frequency sqrt(fx^2 + fy^2)
% Lr : cut-off distance
function S = spectral_density_bandpass_2d(fr,fc,p,normalize,q)
	if (nargin()<4)
		normalize = 1;
	end
	fc = fc/(2*2^(1/3) - 1)^(1/2);
	if (nargin()<5)
		q = 1;

	switch (q)
	case {1}
		% from recombining a first order lowpass filter
		S = (4*(1 - (1 + (fr./fc).^2).^(3/4))).^2./(1 + (fr./fc).^2).^3;
	case {2}
		% from recombining a second order lowpass filter
		S = (4*(1 - ((fr./fc).^2 + 1).^(3/2)).^2)./(1 + (fr/fc).^2).^6;
	otherwise
		error('not yet implemented');
	end

	if (nargin()>2 && ~isempty(p))
		S = S.^p;
	end
	switch(normalize)
		case {-1,0} % nothing
		case {1}
			Sc = spectral_density_bandpass_2d_scale(fc,p);
			S = Sc.*S;
	end % switch
end

