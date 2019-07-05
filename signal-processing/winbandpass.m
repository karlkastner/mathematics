% Mi 24. Feb 09:40:55 CET 2016
% Karl Kastner, Berlin
%% filter with bandpass
function [win t] = winbandpass(dt,T)
	mode = 2;
	switch (mode)
	case {1}
	T = 3/2*T;
	% TODO odd even correction
	t = (0:dt:T)';
%	clear y; L=3; n=1e3; x = linspace(0,3,n)'; y = [sin(pi*x*3/L), sin(pi*x/L) ]; y(:,3) = y(:,1).*y(:,2); plot(x,y); sum(y)
	p = 0;
	win =        sin(pi*t*3/T).*sin(pi*t/T).^p ...
	      + 1.0i*cos(pi*t*3/T).*sin(pi*t/T).^p;
	% TODO, this is not exact
	win = win/norm(win);
%	win = 1/1.584*win;
	case {2}
		p = 1;
		T2 = 1.5*T;
		t = (0:dt:T2-dt)';
		win =       sin(2*pi*t/T).*sin(pi*t/(T2)).^p ...
		      +1.0i*cos(2*pi*t/T).*sin(pi*t/(T2)).^p;
		win = win/norm(win);
	case {3}
		T2 = T;
		t = (0:dt:T2-dt)';
		win =       1*sin(2*pi*t/T) ...
		      + 1i*cos(2*pi*t/T);
		win = win/norm(win);
	end
end

