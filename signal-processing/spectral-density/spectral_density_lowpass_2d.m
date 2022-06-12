% Fri 22 Apr 14:05:05 CEST 2022
function S = sd_lowpass_2d(fx,L,p)
	fx = 2*pi*fx;
	S = (L.^2*fx.^2 + 1).^(-3/2);
	S = S.^(p);
end
	
