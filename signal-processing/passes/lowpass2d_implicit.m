% Thu 24 Jun 16:27:03 CEST 2021
%
% lowpass filter the surface x by solving the implicit relation
% note : this is computationally inefficient and serves demonstration
%
%
function [y] = lowpass2d_implicit(x,rho,a,order,direct)
	if (any(rho >= 1))
		warning('rho must be smaller 1');
	end
	if (length(rho) < 2)
		rho = [rho,rho];
	end
	if (nargin()<4)
		order = 1;
	end
	if (nargin()<5||isempty(direct))
		direct = false;
	end
	n  = size(x);
	tol = [];
	maxit = max(n);
	I = speye(prod(n));
	% y[i,j] = rho*(y[i-1] + y[i+1] + y[j-1] + y[j+1]) + (1-4*rho)*x[i] 
	% (1 - 4 rho)y = rho*[1,..,1,-4,1,...,1]*y + (1-4 rho)*rho*x[i]
	% (1 - 4 rho)y = rho*D2*y + (1-4 rho)*x[i], D2 = D2x+D2y
	% ((1 - 4*rho)*I - rho*D2*y = (1-4*rho)*x
	% y = (1-4*rho)*(((1-4*rho)*I - rho*D2)^-1*x
	% y = ((I - rho/(1-4*rho)*D2)^-1*x
	% lowpass
	%rho = rho./(1-2*(rhox+rhoy)+rhox*rhoy);
	%r = rho./(1 - 2*rho + rho.*rho);
	%r = rho; %./(1 - 2*rho + 2*rho.*rho)
	% TODO this approximation is only exact in the limit rho -> 1
	% Note that the difference approach does not yield r^|-k| in 2D
	% this neither seems to improve with higher order finite differences
	% TODO, try with spectral accurracy
	% or larger domains
	% s = sum(sum(rho^sqrt(i^2+j^2)))
	k = (-200:200);
	s = 2*pi/log(rho(1)).^2
	s = sum(sum(rho(1).^hypot(k,k')))
	r = 1/4*(s-1)./(1-rho(1));
	r = [r,r];

%	switch (direct)
%	case {0,1}

	% order of accuracy for finite difference approximation
	norder = 2;
	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,n-1,norder,'circular');         

	rD2  = (r(1)*D2x + r(2)*D2y);

	if (nargin()>2 && ~isempty(a))
		s = sin(a);
		c = cos(a);
		D2 = c*c*D2x + 2*s*c*Dxy + s*s*D2y;
		% TODO asymmetric scaling
		rD2 = r(1)*D2;
	end

	A = (I - rD2);

%	case {'fourier'}
		[Dx,Dy,D2x_f,Dxy,D2y_f] = fourier_derivative_matrix_2d(n,n-1);
%	end % switch

	y = x;
	for idx=1:order
	switch (direct)
	case {true}
		y = A \ flat(y);
	case {false}
		y = pcg(A,flat(y),tol,maxit);
	case {'fourier'}
		y0 = pcg(A,flat(y),tol,maxit);
		%y0 = y;
		y0 = 1e-3*randn(size(y0))+y0;
		y = pcg(@Afun_fourier,flat(y),tol,maxit,[],[],flat(y0)); %,tol,maxit);
	end % switch
	end % for

	y = reshape(y,n);

	function x = Afun_fourier(y)
		% TODO is this possible with fft2 or are 2 independent fft1 requires?
		%x = y + r*flat(ifft2(D2x_f*D2y_f*fft2(reshape(y,n))));
if (0)
		y_ = fft2(reshape(y,n));
		y_ = flat(y_);
		dy = D2x_f*D2y_f*y_;
		dy = reshape(dy,n);
		dy = ifft2(dy);
		dy = flat(dy);
else
		y_ = fft(reshape(y,n));
		y_ = flat(y_);
		dy = D2y_f*y_;
		dy = reshape(dy,n);
		dy = ifft(dy);
		% second dim
		dy = fft(dy.').';
		dy = flat(dy);
		dy = D2x_f*dy;
		dy = reshape(dy,n);
		dy = ifft(dy.').';
		dy = flat(dy);
end
		x = y - r(1)*dy;
	end
end

