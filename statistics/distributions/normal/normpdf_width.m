function w=normpdf_width(mu,s,p)
	if (nargin()<p)
		p = 0.5;
	end
	w = sqrt(-8*log(p))*s;
end

