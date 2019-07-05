% Sun 18 Dec 21:13:46 CET 2016
% Karl Kastner, Berlin
%
%% binary search on a line
function c=binsearch(x,x0,l,r)
	while (l ~= r)
		c = floor((l+r)/2);
		if (x(c)<x0)
			r = c;
		else
			l = c;
		end
	end
end

