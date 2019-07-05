% Fri 10 Nov 16:21:57 CET 2017
%% van Leer limiter
function phi = vanLeer(theta)
	phi = (theta + abs(theta))./(1+abs(theta));
end

