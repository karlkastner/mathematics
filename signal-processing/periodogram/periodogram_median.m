% Tue 21 Sep 19:37:17 CEST 2021
function S = median_periodogram(Shat,m)
	Shat = fftshift(Shat);
	S  = medfilt1(Shat,m) / (1/2*chi2inv(0.5,2));
	S  = S*1.05;
	S  = ifftshift(S);
end

