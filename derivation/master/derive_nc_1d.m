% Sun Jul 15 12:31:35 MSK 2012
% Karl Kästner, Berlin

% derive closed newton cotes formulae
function [w q] = derive_nc_1d(n)

% quadrature point coordinates
q = (0:n-1)'/(n-1)

% quadrature point vandermonde matrix
A = vander_1d(q,n)';
% coefficients of the polynomial test functions
C = inv(vpa(A))
% coefficients of the integrated polynomial test functions
F = integrate_1d(C)
% evaluate integrated test functions for the constant 1 function
% gives weights for the quadrature points
% w = F(1) - F(0) = F(1)
w = vpa((F*ones(n+1,1))')
% contrary to the irrational Gauss weights and coordinates,
% are the N-C quadrature weights points are exact rational fractions
%w = double(w)
%[num den] = rat(w)
%l=den(1);
%for idx=2:length(den)
%	l = lcm(l,den(idx));
%end	
%l
%w*l
%test
%norm(w-num./den)
whos

%syms x, f1 f2 f3 f4 f5 f6
%F = [f1 f2 f3 f4 f5 f6];
%X = [1 x x^2 x^3 x^4 x^5];
%f = Fv(1:n)*C*X.'
%F = int(f,x)
%I = subs(F,x,1) - subs(F,x,0)

end

function V = vander_1d(x,n)
	l = size(x,1);
	V = zeros(l,n);
	if (n > 0)
		V(:,1) = ones(l,1);
	end
	if (n > 1)
		V(:,2) = x;
	end
	for idx=3:n
		V(:,idx) = V(:,idx-1).*x;
	end
end

% integrates once an 1D polynomial of power n
% rows of C contain coefficients for [1 x x^2 ...]
function I = integrate_1d(C)
	n1 = size(C,1);
	n2 = size(C,2);
	I = vpa(zeros(n1,n2+1));
	I(:,2:n2+1) = (ones(n1,1)*(1./(1:n2))).*C
end % integrate 1d

