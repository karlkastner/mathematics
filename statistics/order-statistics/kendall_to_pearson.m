% Sa 1. Aug 13:16:03 CEST 2015
% Karl Kastner, Berlin
%
%% convert kendall rank correlation coefficient to the person product moment
%% correlation coefficient
%%
%% c.f. Kruska, 1985
function rho = kendall_to_pearson(tau)
	rho = sin(pi*tau/2);
end

