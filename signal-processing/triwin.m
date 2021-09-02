% Fri Mar  6 18:13:35 CET 2015
% Karl Kastner, Berlin
%% triangular filter window
% function w = triwin(x,x0,L)
function w = triwin(x,x0,L)
	if (nargin()<2)
		x0 = 0.5*(x(1)+x(end));
	end
	if (nargin() < 3)
		L = x(end)-x(1);
		% scale to make weight at first and last point non-zero
		nx = length(x);
		L = L*(nx)/(nx-1);
	end
	x = 0.5-abs(x-x0)/L;
	w = x.*(x>0);
	w = w/sum(w);
end

