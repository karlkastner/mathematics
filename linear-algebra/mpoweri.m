% Fri 14 Jul 16:14:25 CEST 2017
%% approximation of A^p, where p is not integer by quadtratic interpolation
function x = mpoweri(A,p,x)
	if (p<0)
		p
		error('here');
	end

	pi_ = floor(p);
	r  = p-pi_;

	x0 = x;
	for idx=1:pi_
		x0 = A*x0;
	end
%	x0 = A^pi*x;
	x1 = A*x0;
	x2 = A*x1;
	x3 = A*x2;

% quadratic	
%	x = (r^2/2 - (3*r)/2 + 1)*x0 + (- r^2 + 2*r)*x1 + (r^2/2 - r/2)*x2;
% cubic: note, this is still problematic, derivatives are discontinuous,
% 	 use better hermite polynomials
	w = [ (- r^3/6 + r^2 - (11*r)/6 + 1); ...
	      (  r^3/2 - (5*r^2)/2 + 3*r); ...
	      (- r^3/2 + 2*r^2 - (3*r)/2); ...
	      (  r^3/6 - r^2/2 + r/3) ];
	x = [x0, x1, x2, x3]*w;
end

