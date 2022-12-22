% Wed 16 Nov 10:19:09 CET 2022
%
% restore the truncated negative half plane of a Fourier series by mirroring
% (F(-k) = conj(F(k)))
%
% the number of tapps of the restored series can be even or odd, default is even
%
% function spectral_density_complete_negative_half_plane(k,f,makeeven)
function [k,f] = fourier_complete_negative_half_plane(k,f,makeeven)
	if (nargin()<3)
		makeeven = true;
	end
	if (k(end)<0)
		% negative half plane already there, nothing to do
	else
		k = cvec(k);
		if (isvector(f))
			f = cvec(f);
		end
		if (makeeven)
			% even
			k = [k; -2*k(end)+k(end-1); -flipud(k(2:end-1))];
			f = [f; 0; conj(flipud(f(2:end-1,:)))];
		else % odd
			% why -1?
			k = [k; -flipud(k(2:end-1))];
			f = [f; conj(flipud(f(2:end-1,:)))];
		end
	end
end % fourier_complete_negative_half_plane

