% Tue Jan 14 12:56:55 WIB 2014
function isnight(t)
	night =   mod(t,1) < 6/24 ...
		| mod(t,1) > 18/24;
end

