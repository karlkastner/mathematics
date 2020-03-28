% 2016-03-08 15:01:45.823753709 +0100
function H = cdfS(h)
	 H = cumsum(h,2) - 0.5*h;
end

