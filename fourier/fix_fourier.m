% Thu  4 May 13:03:17 CEST 2017
% Karl Kastner
%
%% fill gaps (missing data) by means of fourier extrapolation
%%
%% fix periodic data series with fourier interpolation
%% longest gap should not exceed 1/2 of the shortest time span of interest (1/cutoff frequency)
%%  note: this limit equals the position of first side lobe of the ft of a rectangular window with gap length
function [x, Fi, c] = fix_fourier(x,m)
	m   = ceil(m);
	fdx = isfinite(x);
	% compute the truncated fourier transform
	[f, Fi] = nanfft(x,m);
	% expand the missing values
	mdx = ~fdx;
	ir  = isreal(x);
	x(mdx) = Fi(mdx,:)*f;
	if (ir)
		x = real(x);
	end
	if (nargout() > 2)
		c = cond(Fi(fdx,:));
	end
end

