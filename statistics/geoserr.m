% Thu 26 Mar 12:47:57 +08 2020
function serr = geoserr(x)
	n = length(x);
	serr = geostd(x)./sqrt(n);
end

