% Thu 26 Mar 12:49:04 +08 2020
function std_ = geostd(x)
	n     = length(x);
	std_  = geomean(x).*std(log(x));
end
