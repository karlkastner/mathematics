% 2021-10-04 10:49:19.880757071 +0200
function [phi,sphi] = average_angle(a,b)
	if (iscomplex(a))
		b = imag(a);
		a = real(a);
	else if (nargin()<2)
		b = cos(a);
		a = sin(a);	
	end

	sa = std(a);
	sb = std(b);
	
	a = mean(a);
	b = mean(b);

	phi = atan2(b,a);
	sphi = sqrt( (a^2*sb^2 + b^2*sa^2)./(a^2+b^2) );
end
