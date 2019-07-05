% Sun 10 Jul 14:21:16 CEST 2016
% Karl Kastner, Berlin
%% kaiser filter window
function win = kaiserwin(x,x0,L,alpha)
	if (nargin() < 4 || isempty(alpha))
		alpha = 3;
	end
	if (nargin() < 2 || isempty(x0))
		x0 = 0.5*(x(end)+x(1));
	end
	x = x-x0;
	if (nargin() < 3 || isempty(L))
		L = x(end);
	end
	x = x/L;
	% x = (0:N-1)/(N-1);
	%win = besseli(0,pi*alpha*sqrt(1-(2*x-1)^2))./besseli(0,pi*alpha);
	win = besseli(0,pi*alpha*sqrt(1-x.^2))./besseli(0,pi*alpha);
	win(x<-1 | x>1) = 0;
	win = win/sum(win);
end

