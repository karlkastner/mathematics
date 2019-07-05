% 2018-02-02 16:02:59.891927504 +0100
%
%% iir-low pass
function x = lp_iir(rho,x)
	s  = (1-rho);
	x = cvec(x);
	x = filter(s,[1 -rho],x);
end

