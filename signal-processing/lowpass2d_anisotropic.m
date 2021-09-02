% Sun 11 Jul 16:03:45 CEST 2021
% note : this is only well approximated, when the matrices D2
% are approximated by spectral operators
function y = lowpass2d_anisotropic(x,rho,a,mode, order);
	n = size(x);
	nn = prod(n);

if (1)
%	[Dx_l,Dy_l,D2x,Dxy,D2y] = derivative_matrix_2d(n*[1,1],(n-1)*[1,1],-1,'circular');
%	[Dx_r,Dy_r,D2x,Dxy,D2y] = derivative_matrix_2d(n*[1,1],(n-1)*[1,1],1,'circular');
	[Dx_c,Dy_c,D2x_c,Dxy_c,D2y_c] = derivative_matrix_2d(n,(n-1),order,'circular');
%	[Dx_2,Dy_2,D2x_2,Dxy_2,D2y_2] = derivative_matrix_2d(2*n*[1,1],2*(n-1)*[1,1],-1,'circular');
	I=speye(nn);

	c = cos(a);
	s = sin(a);

	% first axis 
	D2s = c*c*D2x_c + c*s*(Dx_c*Dy_c+Dy_c*Dx_c) + s*s*D2y_c;
	% second axis
	D2n = s*s*D2x_c - c*s*(Dx_c*Dy_c+Dy_c*Dx_c) + c*c*D2y_c;

	% renormalize
	rho = rho./(1 - 2*rho + rho.*rho)

	% low-pass matrix
	A  = (I - rho(1)*D2s - rho(2)*D2n);
else
%	D2a  = diffusion_matrix_2d_anisotropic(n,n-1,rho*[1,0],-deg2rad(a(idx)));
	D2  = diffusion_matrix_2d_anisotropic(n,n-1,rho,-a);
	A   = (I - D2a_);
end

	y = flat(x);

	switch (mode)
	case {'low'}
		y = A\y;
	case {'high'}
		y = y - A\y;
	case {'band'}
		y = A\y;
		y = y - A\y;
	otherwise 
		error('')
	end

	y   = reshape(y,n);
%	y_2 = reshape(y_2,2*n,2*n);
%	y_2 = imrotate(y_2,a(idx),'nearest','crop');
%	y_2 = y_2(n/2:end-n/2,n/2:end-n/2);
end


