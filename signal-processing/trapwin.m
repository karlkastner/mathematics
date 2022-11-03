% Fri 17 Feb 18:55:33 CET 2017
%% trapezoidal filter window
% p = 0 : rectangular
% p = 1 : triangular
% function w = trapwin(x,x0,L,p)
function w = trapwin(x,x0,L,p)
	if (nargin()<2 || isempty(x0))
		x0 = 0.5*(x(1)+x(end));
	end
	if (nargin() < 3 || isempty(L))
		L = x(end)-x(1);
		% scale to make weight at first and last point non-zero
		nx = length(x);
		L = L*(nx)/(nx-1);
	end
	if (nargin() < 4)
		p = 0.5;
	end
	w = zeros(size(x));
	x = (x-x0)/L;
	fdx =  x >-0.5 & x<-0.5+0.5*p;
	w(fdx) = 2*(x(fdx)+0.5)/(p);
	fdx = x >= -0.5 + 0.5*p & x <= 0.5-0.5*p;
	w(fdx) = 1;
	fdx = x > 0.5-0.5*p & x < 0.5;
	w(fdx) = -2*(x(fdx)-0.5)/(p);

%	a=asin((2*tukeywin(100,0.5))-1)/pi+0.5
%	trapmf(x,[0 25 75 100])

%	x = 0.5-abs(x-x0)/L;
%	w = x.*(x>0);
%	w = w/sum(w);
end

