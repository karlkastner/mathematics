% Sun 18 Sep 10:36:54 CEST 2016
%% danielle fourier window
function w = danielle_window(n)
	a =-(2^(-1/2)*(n^2 - n)^(1/2) - n + 1)/(n - 2);
	% note: this is asymptoically 1-1/sqrt(2) ~ 0.29
	w = [a ones(1,n-1) a];
end

