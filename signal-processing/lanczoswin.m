% Thu  9 Feb 11:58:31 CET 2017
% Karl Kastner, Berlin
%% Lanczos window
function win = lanczoswin(x,x0,L,a)
	if (nargin()<2 || isempty(x0))
		x0 = 0.5*(x(end)+x(1));
	end
	if (nargin()<3 || isempty(L))
		L = 0.5*(x(end)-x(1));
	end
	if (nargin()<4)
		a = 2;
	end
	% shift and normalise
	x = (x-x0)/L;

	if (0)
		A   = vander_1d(cvec(x),[1 0 1]);
		w   = hanwin(x);
		%A(:,2) = [];
		A   = diag(w)*A;
		W   = (A'*A)\(A');
		win = W(1,:)';
	else
		win = a*sin(2*pi*x).*sin(2*pi*x/a)./(4*pi^2*x.^2);
		if (mod(length(x),2)==1)
			win((length(x)+1)/2) = 1;
		end
		win(x<-1) = 0;
		win(x>+1) = 0;
	end
	win = win/sum(win);
end % lanczoswin

