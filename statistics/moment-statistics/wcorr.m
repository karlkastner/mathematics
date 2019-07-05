% Thu 10 Aug 18:27:54 CEST 2017
%
%% correlation of two vectors when samples are weighted
function [rho_12, v_12, v_1, v_2] = wcorr(w1,x1,w2,x2)
	v_12 = wcov(w1,x1,w2,x2);
	v_1  = wvar(w1,x1);
	v_2  = wvar(w2,x2);
	rho_12 = v_12./sqrt(v_1.*v_2);
end

