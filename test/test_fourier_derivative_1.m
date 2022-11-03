 L=1;
 p=1;
 n=1e3;
 x = innerspace(0,L,n)';
 o = 2*pi;
 y=sin(o*x);

D = derivative_matrix_1_1d(n,L,2,'circular');
Di = D;
clf
for idx=1:4
	subplot(2,2,idx)
	dy_dx      = ifft(fourier_derivative_matrix_1d(n,L,idx)*fft(y));
	switch (idx)
	case {1}	
		dy_dx(:,2) = o*cos(2*pi*x);
	case {2}	
		dy_dx(:,2) = -(o)^2*sin(2*pi*x);
	case {3}
		dy_dx(:,2) = -o^3*cos(2*pi*x);
	case {4}
		dy_dx(:,2) = o^4*sin(2*pi*x);
	end

	dy_dx(:,3) = Di*y;
	plot(x,[y,dy_dx]);
%	r  = median(dy_dx(:,1)./dy_dx(:,2))
	Di = D*Di;
end
