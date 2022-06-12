% Tue 30 Nov 16:57:02 CET 2021
% note : this is the one-sided filter,
%        the equation is only valid for p = 1
function S = spectral_density_lowpass_one_sided(j,r,p,n)
	dk = 2*pi/n;
	S = (  1 - r*cos(dk*j) - r^n*(1 - r*cos(dk*j))) ...                                
            ./((1 - 2*r*cos(dk*j) + r^2));
	S = S.^p;
end
