% Sat Aug 11 02:08:59 MSK 2012
% Karl Kästner, Berlin

% last nonzero partial derivatives of the polynomial
function D = derivative_2d(C, p)

	f = zeros(p+1,1);
	f(1) = 1;
	for pdx=1:p
		f(pdx+1) = pdx*f(pdx);
	end
	
	n0 = p*(p+1)/2;
	for idx=1:size(C,1)
		ndx=1;
		for dx=p:-1:0
			dy=p-dx;
			D(idx,ndx) = f(dx+1)*f(dy+1)*C(idx,n0+ndx);
			ndx = ndx+1;
		end
	end
end % derivative_2d

