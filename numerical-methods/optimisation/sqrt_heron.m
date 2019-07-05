% Mon 18 Jul 18:07:54 CEST 2016
%
%% Heron's method for the square root
function xi = heron(x,n)
	if (nargin < 2)
		n = 10;
	end
	xi = x;
	S  = x;
	for idx=1:n
		xi = 1/2*(xi + S/xi);
	end	
end

