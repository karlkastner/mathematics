% Mon 12 Mar 15:18:59 CET 2018
% Mon 12 Mar 14:10:46 CET 2018
%
%% derivative matrix on a curvilinear grid
%
% for higher order see KWOK 1985
% Wiesenberger 2017 for metric tensor G elements
% BERNARD 1991 for simplest explanation of derivatives
% Thompson 1980 for simplest explanation of second derivatives
% awa Tensors and the Equations of Fluid Motion (chap5_3.pdf)

% on a circle (polar):
%grad phi e_theta = 1/r dphi/dtheta
%grad phi e_r = dphi/dr
function [D1x, D1y, D2x, Dxy, D2y, L] = derivative_curvilinear(x,y,isorthogonal)
	if (nargin()<3)
		isorthogonal = true;
	end

	Dxy = [];

	n = size(x);

	% unit difference matrices
	[Ds,Dn,D2s,Dsn,D2n] = derivative_matrix_2d(n,n-1);

	% forward transformation matrix and step width
	dx_ds = Ds*flat(x);
	dx_dn = Dn*flat(x);
	dy_ds = Ds*flat(y);
	dy_dn = Dn*flat(y);

	d2x_ds2 = D2s*flat(x);
	d2x_dsn = Dsn*flat(x);
	d2x_dn2 = D2n*flat(x);
	d2y_ds2 = D2s*flat(y);
	d2y_dsn = Dsn*flat(y);
	d2y_dn2 = D2n*flat(y);

	J(1,1,:) = dx_ds;
	J(1,2,:) = dy_ds;
	J(2,1,:) = dx_dn;
	J(2,2,:) = dy_dn;
	% inverse transformation matrix
	if (1)
		Ji = inv2x2(J);
	else
		Ji = J;
		J  = inv2x2(Ji);
	end
	% [ ds/dx, dn/dx] du/ds = [ ] du/dx
	% [		] du/dn   [ ] du/dy
	ds_dx = squeeze(Ji(1,1,:));
	dn_dx = squeeze(Ji(1,2,:));
	ds_dy = squeeze(Ji(2,1,:));
	dn_dy = squeeze(Ji(2,2,:));


	J11   = squeeze(J(1,1,:));
	J12   = squeeze(J(1,2,:));
	J22   = squeeze(J(2,2,:));
	J21   = squeeze(J(2,1,:));
	Ji11  = squeeze(Ji(1,1,:));
	Ji12  = squeeze(Ji(1,2,:));
	Ji21  = squeeze(Ji(2,1,:));
	Ji22  = squeeze(Ji(2,2,:));
	sdetJ = sqrt(J11.*J22);

	

	% d^2/dx^2 f(s,n) =    ds^2/dx^2 df/ds 
        %                   + (ds/dx)^2  d^2f/ds^2 
        %                   + (dn/dx)^2  d^2f/dn^2                         
        %                   +  dn^2/dx^2 d^2f/dn^2                         
	%                   + diag(sparse(dn_dx))*Dn;

	% page 244 in Numerical Solution of the Incompressible Navier-Stokes Equations, Quartapelle
	% r : position of point on curve
	% hs = dr/ds
	% A-6.2 in 3A978-94-009-8352-6%2F1
	hs = flat(1./abs(hypot(cdiff(x),cdiff(y))));
	% hn = dr/dn
	hn = flat(1./abs(hypot(cdiff(x'),cdiff(y'))'));

	% in orthogonal coordinates : g^ii = 1/g_ii, g^ij = 0, i neq j

	% in orthogonal curvilinear coordinates the off-diagonal elements vanish (p20,ch 15)
	% ch15p11
	%  grad f = sum_ij g^ij df/du ei
	%  
	% simplifies in curvilinear coordinates to
	% 

	if (1) % ~isorthogonal)
		
		% c.f. A-5.1 3A978-94-009-8352-6
		% without orthongonality ch15p11
			% D1x = diag(sparse(dy_dn./j))*Ds + diag(sparse(-dy_ds./j))*Dn;
			% D1y = diag(sparse(dx_ds./j))*Dn + diag(sparse(-dx_dn./j))*Ds;
		% from matrix inverse identity
		% swap a11 and a22, invert sign of a12 and a21 :
		% dJ   = dx_ds.*dy_dn-dx_dn.*dy_ds;
		% ds_dx :=  dy_dn/dJ = Ji11
		% dn_dx := -dy_ds/dJ = Ji12
		% ds_dy :=  dx_ds/dJ = Ji21
		% dn_dy := -dx_dn/dJ = Ji22
		% hs = sqrt(dx_ds^2+dy_ds)^2
		% d/dx f(s,n) = ds/dx df/ds  + dn/dx df/dn
		D1x = diag(sparse(ds_dx))*Ds + diag(sparse(dn_dx))*Dn;
		D1y = diag(sparse(ds_dy))*Ds + diag(sparse(dn_dy))*Dn;

		if (0)
			% A-6.3, this seems not correct, as the working eq above,
			%        dx is formed by dy
			hs = squeeze(Ji(1,1,:));
			hn = squeeze(Ji(2,2,:));
			D1x =    diag(sparse(hs.^2.*dx_ds))*Ds ...
		      	       + diag(sparse(hn.^2.*dx_dn))*Dn;
			D1y =   diag(sparse(hs.^2*dy_ds))*Ds ...
		              + diag(sparse(hn.^2*dy_dn))*Dn;
		end
	%	D1x = diag(sparse(1./(hs.*hn)))*...
	%		(   Dn*diag(sparse(hs./hn))*Dn ...
	%		  + Ds*diag(sparse(hs./hn))*Ds );

	else
		% with orthogonal mesh identity
		% ch15p20
		D1x = diag(sparse(dx_ds))*Ds;	
		D1y = diag(sparse(dy_dn))*Dn;	
	end


	if (0) %~isorthogonal)

	else

	d2s_dx2 = D1x*ds_dx;
	d2n_dx2 = D1x*dn_dx;
	d2s_dy2 = D1y*ds_dy;
	d2n_dy2 = D1y*dn_dy;

	D2x =   diag(sparse(d2s_dx2))*Ds ...
	      + diag(sparse(ds_dx.^2))*D2s ...
	      + diag(sparse(dn_dx.^2))*D2n ...
	      + diag(sparse(d2n_dx2))*Dn;
	D2y =   diag(sparse(d2s_dy2))*Ds ...
	      + diag(sparse(ds_dy.^2))*D2s ...
	      + diag(sparse(dn_dy.^2))*D2n ...
	      + diag(sparse(d2s_dy2))*Dn;
	L = D2x+D2y;
	% L = D1x*D1x + D1y*D1y;

if (0)
dx_ds.^2+dy_ds.^2
dy_ds.*dy+dy_ds.^2
'honk'
pause

	L = (   (ds_dx.^2 + ds_dy.^2)*D2s ...
	      + 2*(ds_dx*dn_dx + ds_dy*dn_dy)*DsDn ...
	      +   (dn_dx.^2 + dn_dy.^2)*D2n ...
	      +   (a)*Ds ...
	      +   (a)*Dn ...
	    );
end

	if (1)		
		% ch15p23, non-orthogonal grid
		J = sqrt(J11.*J22);
		L_ = diag(sparse(1./J)) ...
		* ( ...
			  Ds*diag(sparse(J.*Ji11))*Ds ...
			+ Ds*diag(sparse(J.*Ji12))*Dn ... 
			+ Dn*diag(sparse(J.*Ji21))*Ds ...
			+ Dn*diag(sparse(J.*Ji22))*Dn ...
		 );

		x_n = diag(sparse(dx_dn));
		x_s = diag(sparse(dx_ds));
		x_ss = diag(sparse(d2x_ds2));
		x_sn = diag(sparse(d2x_dsn));
		x_nn = diag(sparse(d2x_dn2));

		y_n = diag(sparse(dy_dn));
		y_s = diag(sparse(dy_ds));
		y_ss = diag(sparse(d2y_ds2));
		y_sn = diag(sparse(d2y_dsn));
		y_nn = diag(sparse(d2y_dn2));

		Ji = 1./(dx_ds.*dy_dn - dx_dn.*dy_ds);
		Ji = diag(sparse(Ji));

		L = 1*Ji^2*(x_n^2 + y_n^2)*D2s  ...
		  - 2*Ji^2*(x_n*x_s + y_n*y_s)*Dsn  ...
		  + Ji^2*(x_s^2 + y_s^2)*D2n  ...
		 - (  (      x_s*(x_nn*Ji + (x_n*(x_n*y_sn + x_nn*y_s - x_s*y_nn - x_sn*y_n))*Ji^2) ... 
			 -   x_n*(x_sn*Ji + (x_n*(x_n*y_ss - x_ss*y_n - x_s*y_sn + x_sn*y_s))*Ji^2))*Ji .... 
		    + (      y_s*(y_nn*Ji + (y_n*(x_n*y_sn + x_nn*y_s - x_s*y_nn - x_sn*y_n))*Ji^2) ...
			 - y_n*(y_sn*Ji + (y_n*(x_n*y_ss - x_ss*y_n - x_s*y_sn + x_sn*y_s))*Ji^2))*Ji)*Ds ...
		+  (  (    x_s*(x_sn*Ji + (x_s*(x_n*y_sn + x_nn*y_s - x_s*y_nn - x_sn*y_n))*Ji^2) ...
			 - x_n*(x_ss*Ji + (x_s*(x_n*y_ss - x_ss*y_n - x_s*y_sn + x_sn*y_s))*Ji^2))*Ji ... 
		    + (    y_s*(y_sn*Ji + (y_s*(x_n*y_sn + x_nn*y_s - x_s*y_nn - x_sn*y_n))*Ji^2) ...
			 - y_n*(y_ss*Ji + (y_s*(x_n*y_ss - x_ss*y_n - x_s*y_sn + x_sn*y_s))*Ji^2))*Ji)*Dn;

		% ch15p23, orthogonal-grid
		%j     = dx_ds.*dy_dn-dx_dn.*dy_ds;
		%sdetJ = 1./j;
		if (0)
		sdetJ = sqrt(squeeze(J(1,1,:)).*squeeze(J(2,2,:)));
		L = diag(sparse(1./sdetJ)) ...
		    * ( Ds*diag(sparse(sdetJ./squeeze(J(1,1,:))))*Ds ...
		      + Dn*diag(sparse(sdetJ./squeeze(J(2,2,:))))*Dn );
		end
	end

end % derivative_curvilinear

