% Fri 24 Jan 19:25:30 +08 2020
	syms x1 x2 x3 x4
	clear f g
	x = [x1;x2;x3;x4];
	f = [     x(1);
                  (x(2)-1).^3;
                  x(3).^2;
		  (x(4)+0.5).^4
	    ]
	for idx=1:length(x)
		g(idx,1) = diff(f(idx),x(idx));
	end
	f = matlabFunction(f,'var',x)
	g = matlabFunction(g,'var',x)
	
	f = @(x) f(x(1),x(2),x(3),x(4));
	g = @(x) g(x(1),x(2),x(3),x(4));
	
	x = rand(size(x));
	x0 = fzero_newton(f,g,x)

	x0(:,2) = fzero_bisect(f,-2*ones(size(x)),1*ones(size(x))+rand(size(x)))

