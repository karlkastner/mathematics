% Wed 18 Jan 09:51:28 CET 2017
function [f,n] = window_make_odd(f,n)
	f = rvec(f);
	if (0 == mod(length(f),2))
	if (nargin()<2)
	n = 1;
	end
	switch (n)
	case {1} % linear
		f = 0.5*([0,f] + [f,0]);
	case {3}
		f = 1/48*(-3*[f,0,0,0] + 27*[0,f,0,0] + 27*[0,0,f,0] -3*[0,0,0,f]);
	case {5}
		f = 1/256*(    3*[f,0,0,0,0,0] ...
			   -  25*[0,f,0,0,0,0] ...
			   + 150*[0,0,f,0,0,0] ...
			   + 150*[0,0,0,f,0,0] ...
			   -  25*[0,0,0,0,f,0] ...
			   +   3*[0,0,0,0,0,f]);
	end
end

