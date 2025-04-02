

% S = s*1/f^a
% R = s*2 int_0^infty 1/f^a exp(-2 i pi f x) df 
%   ~     1/x^(1-a)

function acf = poweracf(x,a,dim)
	if (nargin()<3)
		dim = 1;	
	end
	scale = (2*pi)^(a-1);
	acf   = scale./(x.^(1+a));
end

