function x = make_monotonic(x,sig)
	for idx=2:length(x)
		if (    ((sig > 0) & (x(idx)<=x(idx-1))) ...
		     || ((sig < 0) & (x(idx)>=x(idx-1))) ...
		   )	
			x(idx) = x(idx-1) + sig*eps*(1 + abs(x(idx-1)));
		end
	end
end
