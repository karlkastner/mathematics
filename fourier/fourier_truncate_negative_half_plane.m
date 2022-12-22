% Wed 16 Nov 10:15:25 CET 2022
%
% truncate the redundant negative half plane of a Fourier series of a real Function
% as it is by definition symmetric (f(k) = conj(f(-k)))
function [k,f] = fourier_truncate_negative_half_plane(k,f)
	if (k(end)>=0)
		% negative half plane already truncated, nothing to do
	else
		n = length(k);
		if (isvector(f))
			f = cvec(f);
		end
		if (0 == mod(n,2))
			% even
			k = k(1:n/2);
			f  = conj(f(1:n/2,:));
		else
			% odd
			k = k(1:(n+1)/2);
			f  = conj(f(1:(n+1)/2,:));
		end
	end % else of k(end)>0
end % fourier_truncate_negative_half_plane

