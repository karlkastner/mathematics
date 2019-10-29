function test_df_2d()
	clf
	for d=0:3
	n=2^14;
	k = 2;
	L0 = 10;
	X = L0*(0:n-1)'/(n);
	%X = 2*pi*(1:n)'/n;
	f  = sin(k*pi*X);
	L = 1;
	f = (X-0.5*L0).^3;
	f = exp(-L^2*(X-L0/2).^2);
	Dx = df_2d([n,1],d);
	[void Dy] = df_2d([1,n],d);
	norm(Dx-Dy')
	switch (mod(d,4))
		case {0}
			df = f;
		case {1}
			df = k^2*cos(k*pi*X);
			df = 3*(X-0.5*L0).^2;
			df = -2*L^2*(X-0.5*L0).*f;
		case {2}
			df = -(k^d)*sin(k*pi*X);
			df = 6*(X-0.5*L0);
			df = ((-2*L^2*(X-0.5*L0)).^2 - 2*L^2).*f;
		case {3}
			df = -(k^d)*cos(k*pi*X);
			df = (-2*L^2*(X-0.5*L0)).^2.*f - 2*L^2*f;
			%df = (((-2*L^2*(X-0.5*L0)).^2 - 2*L^2).*(-2*L^2*(X-0.5*L0)) + 8*L^2*(X-0.5*L0)).*f;
			df = ((8*L^4*(X - 0.5*L0)) - (L^2*2*(X - 0.5*L0).*(4*L^4*(X - 0.5*L0).^2 - 2*L^2))).*f;

	end
	f4 = real(ifft(Dx.*fft(f)))/L0^d;
%	f4 = real(ifft(Dx.*fft(f*(L/(2*pi))^d)));
%	df./f4
	norm(df-f4)
	subplot(2,2,d+1)
	plot(X,[f df f4]);
	end
end
