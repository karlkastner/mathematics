% Sun May 20 18:20:41 MSK 2012
% Karl Kästner, Berlin

% derives the nth-order polynomial once
%
% [1 x y x^2 xy y^2]'*[c00 c10 c01 c20 c11 c02]'
% => dv/dx : [1 x y]'*[c10 2*c20 c11]
%    dv/dy : [1 x y]'*[c20 c11 2*c02]
%
function [C_dx, C_dy, Dx,Dy] = derivative_2d(C,n);
	a    = size(C,1);
	b    = n*(n+1)/2;
	if (~issym(C))
		C_dx = zeros(b,a);
		C_dy = zeros(b,a);
	end
	s = 0;
	m = 1;
	for idx=1:n
		s = s+idx;
		for jdx=1:idx
			C_dx(m,:) = (idx+1-jdx)*C(s+jdx,:);
			C_dy(m,:) =         jdx*C(s+jdx+1,:);
			m = m+1;
		end
	end
end % derivative_2d



