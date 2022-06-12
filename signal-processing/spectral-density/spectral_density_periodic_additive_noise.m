% Wed  8 Sep 15:58:26 CEST 2021
%
% spectral density of a periodic function with additive noise
%
% y(t) = (sum a_i cos(2*pi*fc_i*t) + b_i cos(2*pi*fc_i*t)) + e(t)
% e ~ N(0,s)
% 
% fc : frequencies of components
% a, b : amplitude of cos and sin part
% s : standard deviation of noise
% L : domain length
% n : number of samples
% normalize : true for normalization
function [S_fc_mu,S_fc_sd,S_f_mu,S_f_sd] = spectral_density_additive_noise(fc,a,b,se,L,n,normalize)
	if (normalize)
		df  = 1/L;
		den = (1/2*sum(a.^2+b.^2) + 1/2*se.^2)*df;
	else
		den = 1;
	end
	S_fc_mu = (1/2*(a.^2+b.^2) + 1/n*se.^2)/den;
	S_f_mu = 1/n*se.^2/den;

	S_fc_sd = sqrt(2/n*S_fc_mu/den)*se;
	S_f_sd = S_f_mu;
end

