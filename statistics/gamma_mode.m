% 2023-01-15 17:54:09.272681131 +0100
function [fc,Sc] = gamma_mode(a,b,flag)
	if (a<1)
		fc = 0;
	else
	%	fc = a-1/b;
		fc = (a-1)*b;
	end
	if (nargin()<3 || ~flag)
		Sc = gampdf(fc,a,b);
	else
		Sc = gampdf_man(fc,a,b)
	end
end

