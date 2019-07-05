% Fri  8 Jul 16:41:22 CEST 2016
% Karl Kastner, Berlin
%% the flat top window
function w = flattopwin(x,x0,L)

	if (nargin() < 3 )
		L = (x(end)-x(1));
		%x = x/L;
	end
	if (nargin() < 2)
		x0 = (x(1)+x(end))/2;
		x0 = x0 + L/2;
	end
	x = x-x0;
	x = x/L;
		
	% third order
	A3 = [0.2811, -0.5209, 0.1980];
%	A_ = [1;
%	     -1.93;
%	     +1.29;
%	     -0.388;
%	     +0.028 ];
	A5 = [ +0.21557895;
              -0.41663158;
	      +0.277263158;
              -0.083578947;
              +0.006947368];
	A = A3;
	w = A(1)*ones(size(x));
	for idx=2:length(A)
		w = w + A(idx)*cos((idx-1)*2*pi*x);
	end
%	w = a0  - a1*cos(2*pi*x) ...
%       		+ a2*cos(4*pi*x) ...
%		- a3*cos(6*pi*x) ...
%       	        + a4*cos(8*pi*x);
end

