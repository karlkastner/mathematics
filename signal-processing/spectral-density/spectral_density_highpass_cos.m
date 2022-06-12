% 2022-01-12 15:58:19.660486387 +0100
%
%% consine shaped spectral density of a highpass filter
%
% function [S] = spectral_density_highpass_cos(f,fc,p)
function [S] = spectral_density_highpass_cos(f,fc,p)
	f0 = 0.5*pi*fc/acos(1-2*(0.5.^(1/p)));
	S = (1-1/2*(1+cos(0.5*pi*f/f0))).^p;
	S(abs(f)>2*f0) = 1;
end
