% Mon 19 Dec 17:03:02 CET 2022
% gaussian window in 1D
function [w, dof] = gausswin1(n, nf)
	fx = fourier_axis(1,n);

	% Gaussian window
	% 0.5 = exp(-1/2*(nf/2)^2/s^2)
	% -2*log(0.5) = (nf/2)^2/s^2
	% -8*log(0.5) = nf^2/s^2
	% s^2 = -nf^2/(8*log(0.5))
	% s^2 =  nf^2/(8*log(2))
	% s   =  nf/sqrt(8*log(2))
	s = nf/sqrt(8*log(2));
	w = normpdf(fx,0,s);

	% normalize sum of weigths to 1
	w = w/sum(w(:));

	% degrees of freedom of window
	% when values are uncorrelated
	dof = sum(w(:)).^2./(sum(w(:).^2));
end

