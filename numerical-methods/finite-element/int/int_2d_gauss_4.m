% Thu 14 May 19:19:27 +08 2020
% Karl Kästner, Berlin
% Dunavant
% baricentric coordinates for gauss quadrature
% 6 points, 4th order exact, E = O(h^5)
function [w, b, flag] = int_2d_gauss_3()
	% weights
	w = [
	0.223381589678011
	0.223381589678011
	0.223381589678011
	0.109951743655322 
	0.109951743655322 
	0.109951743655322 
	];
	
	% baricentric coordinates
	b = [
	0.108103018168070 0.445948490915965 0.445948490915965
	0.445948490915965 0.108103018168070 0.445948490915965 
	0.445948490915965 0.445948490915965 0.108103018168070 
	0.816847572980459 0.091576213509771 0.091576213509771
	0.091576213509771 0.816847572980459 0.091576213509771 
	0.091576213509771 0.091576213509771 0.816847572980459
	];
	flag = 0;
end
