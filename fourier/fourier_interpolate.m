% Sat  8 Sep 13:56:11 CEST 2018
% Karl Kastner, Berlin
%
%% interpolate samples y sampled at moments (location) t to locations ti

function [yi,c] = fourier_interpolate(t,y,ti,Tf)
	t0 = t(1);
	c  = fourier_fit(Tf,t0,t,y);
	yi = fourier_predict(Tf,t0,c,ti); 
end

