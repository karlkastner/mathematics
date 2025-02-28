% Sat 26 Oct 10:28:18 CEST 2024
% converts roots to polynomial coefficients
% note that there are infinitely many solutions
% scaled by a constant, in particular -1 (flipped upside down)
function c = roots2coeffs(r)
	% p=roots2coeffs([2,3,4])
	% syms x;
	% expand((x-2)*(x-3)*(x-4))
	% roots(p)   
	c = zeros(1,length(r)+1);
	% c = prod(x-ri)
	%     (x-r1)
	%     (x - r2)(x - r1) =  x2*(x - r1) - r2*(x - r1)
	c(1) =   1;
	c(2) = -r(1);
	for idx=2:length(r)
		c(2:idx+1)  = c(2:idx+1) - r(idx)*(c(1:idx));
	end % for idx
end % roots2coeffs

