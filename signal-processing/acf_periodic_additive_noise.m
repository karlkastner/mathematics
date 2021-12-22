% Mon  6 Sep 17:12:01 CEST 2021
% autocorrelation of a cosine with gaussian noise
% a : std of the cosine (sqrt(2)*amplitude)
% s : std of the white noise
function [acf,a_acf] = autocorrelation_cos_noisy(x,a,s)
	a_acf = a.^2./(a.^2+2*s.^2);
	acf   = a_.*cos(x);
	acf(x==0) = 1;
end

