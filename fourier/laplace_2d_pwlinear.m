% Tue 15 May 15:17:12 CEST 2018
%
%% solution to the Laplacian in two dimensions for a finite rectangular domain
%% with piecewise constant boundary conditions
%% linear system with 4 unknowns per freqency component
%% these are coefficients of s,c,sh,ch
%%	(pu*(s + c) + qu*(s' + c'))*(shu + chu) = ru		% upper bc
%%	(pd*(s + c) + qd*(s' + c'))*(shd + chd) = rd		% lower bc
%%	( (sl + cl)*( pl*(shl + chl) + ql*(shl' + chl')) = rl	% left bc
%%	( (sr + cr)*( pr*(shr + chr) + qr*(shr' + chr')) = rr	% right bc
%%
%%  least squares with piecewise integration
%% [x0,p,q,r] piecewise linear polynomials at the boundaries
%
% function [u,v,Phi,coeff] = laplace_2d_pwlinear(x,y,X0,Y0,kmax)
%
% TODO, renormalize x and y to 2pi
function [k,c,x,y,Phi,u,v] = laplace_2d_pwlinear(bc,kmax,L,n)
	% determine coefficients
	c   = zeros(4,kmax);
	for k=0:kmax
		A   = 0;
		b   = 0;
		% set up 4x4 system of equations
		% up
		X0   = bc(1).X;
		p    = bc(1).p;
		q    = 1-p;
		rhs0 = bc(1).rhs;
		for id=1:length(X0)-1
		end
		% down
		X0   = bc(1).X;
		p    = bc(1).p;
		q    = 1-p;
		rhs0 = bc(1).rhs;
		for id=1:length(X0)-1
		end
		% left
		X0   = bc(1).X;
		p    = bc(1).p;
		q    = 1-p;
		rhs0 = bc(1).rhs;
		for id=1:length(X0)-1
		end
		% right
		X0   = bc(1).X;
		p    = bc(1).p;
		q    = 1-p;
		rhs0 = bc(1).rhs;
		for id=1:length(X0)-1
		end
		% solve for coefficients
		c(:,k) = A \ b;
	end
	if (nargout() > 2)
		x = linspace(0,L(1),n(1))';
		y = linspace(0,L(2),n(2))';
		%nx = length(x);
		%ny = length(y);
		Phi = zeros(n(1),n(2));
		% expand phi
		for k=0:kmax
			px =    bsxfun(@times,sin(x*k),c(1,:)) ...                                             
			      + bsxfun(@times,cos(x*k),c(2,:));
			% py = bsxfun(@times,py,1./sinh(k*w))
			py = 	bsxfun(@times,sinh(y*k),c(3,:));
			      +	bsxfun(@times,cosh(y*k),c(4,:));
			Phi = Phi + (px*py');
		end
	end
	if (nargout() > 3)
		% expand u and v, (u', v')
		% TODO, truely expand, no numerical derivatives
		v =  bsxfun(@times, cdiff(Phi')',cdiff(x));
		u =  bsxfun(@times,-cdiff(Phi), cdiff(y));
	end
end % laplace_2d_pwlinear

