% Sa 1. Aug 13:13:35 CEST 2015
% Karl Kastner, Berlin
%
%% conversion of pearson to kendall correlation coefficient
%% c.f. Kruskal 1958
function rho_k = pearson_to_kendall(rho_p)
	rho_k = (2/pi)*asin(rho_p);
end

