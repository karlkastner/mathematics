% Fri 18 Nov 21:46:55 CET 2022
function R = autocorr_forest(a,x)
	R = exp(-a*abs(x));
end

