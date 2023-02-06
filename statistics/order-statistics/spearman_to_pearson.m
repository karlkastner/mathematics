% Sa 1. Aug 15:21:52 CEST 2015
% Karl Kastner, Berlin
%
%% conversion of spearman rank to person product moment correlation coefficient
% c.f. Handbook of Engineering Hydrology: Modeling, Climate Change, and Variability, p. 520, eq. 25.5
% kruskal 1958, eq 6.4
function rho = spearman_to_pearson(rho)
	rho = 2*sin(pi/6*rho);
end

