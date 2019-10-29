	n = 10;
	n2 = n*2;
%	n=n+1;
%	n2=n2+1;

	[A0 X0] = harmonic_oscillator(n,1);
	[A1 X1] = harmonic_oscillator(n2,1);
	
	[v0 e0] = eigs(A0,[],1,'SM'); 
	[v1 e1] = eigs(A1,[],1,'SM'); 

	L = 10;
	L=1;
	
	h = L/(n+1);
	h2 = L/(n2+1);

	X0 = (h*(1:n)'-0.5*L);
	X1 = (h2*(1:n2)'-0.5*L);

	A0 = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	D2 = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	A1 = 1/h2^2*spdiags(ones(n2,1)*[1 -2 1],-1:1,n2,n2);
	D4 = 1/h^4*spdiags(ones(n,1)*[1 -4 6 -4 1],-2:2,n,n);
	

%	x=zeros(n,1); x((n-1)/2+1)=(n);
%	x = [1:(n-1)/2+1 (n-1)/2:-1:1]'.^2;
%	x = 1/2*X0.^2 + 1/factorial(4)*X0.^4; % + 1/factorial(6)*X0.^6 + 1/factorial(8)*X0.^8;
%	x=zeros(n,1); x((n)/2+1)=(n);
%	x = exp(-X0.^4);
%	x = sin(2*pi*X0);
	x = exp(-X0.^2);
	v0 = A0*x;
%	x=zeros(n2,1); x((n2-1)/2+1)=(n2);
%	x = [1:(n2-1)/2+1 (n2-1)/2:-1:1]'.^2;
%	x = 1/2*X1.^2 + 1/factorial(4)*X1.^4; % + 1/factorial(6)*X1.^6 + 1/factorial(8)*X1.^8;
%	x=zeros(n2,1);
%	x((n2)/2+1)=(n2);
%	x = exp(-X1.^4);
%	x = sin(2*pi*X1);
	x = exp(-X1.^2);
	v1 = A1*x;

	norm(v0)
	norm(v1)
	
	vp=interp1(X1, v1, X0,'linear');
%	vp = v1(1:2:end);
%	x_ = [X0(1:end) X1(2:2:end)];
	x_(1:10,:)
%	vp=vp/norm(vp);
%	v0 = v0/norm(v0);
%	v1 = v1/norm(v1);
%	norm(vp - v0) 
	err = vp - v0;
	err_est = -1/12*h^2*D2*v0;
%	err_est = 1/12*h^2*D2*v0;
	norm(err - err_est)
	err_est2 = -1/12*h^2*D2*D2*v0 + 1/90*h^6*D4*D2*v0;
	%norm(err - err_est)
	subplot(2,2,1)
	plot(X0,v0);
	subplot(2,2,2)
	plot(X0,vp);
	subplot(2,2,3)
	plot(X0(3:end-2),[err(3:end-2)])
	subplot(2,2,4)
	plot(X0(4:end-3),[err_est(4:end-3)]); hold on
	subplot(2,2,4)
	plot(X0(4:end-3),[err_est2(4:end-3)],'r'); hold off


	

	A0 = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	eigs(A0,[],1,'SM') +pi^2
	eigs(A0-1/12*h^2*A0^2,[],1,'SM') +pi^2
	eigs(A0-1/12*h^2*A0^2+1/90*h^4*A0^3,[],1,'SM') +pi^2
