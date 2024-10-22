% 2024-09-27 09:24:27.614761185 +0200
% Karl Kastner, Berlin
%
%% elements of the discretization matrix of a two-dimensional
%% reaction-advection-diffusion system`
function di = rad2d_coeffs2diags(d,ad,afun,L,n,nvar,val)
	dt = 1;
	dx = L./n;
	di = zeros(n(1),n(2),4+nvar,nvar);
	for idx=1:nvar
		di(:,:,1,idx) = -dt*d(1)/dx(1)^2;
		di(:,:,2,idx) = -dt*d(1)/dx(1)^2;
		di(:,:,3,idx) = -dt*d(2)/dx(2)^2;
		di(:,:,4,idx) = -dt*d(2)/dx(2)^2;
		for jdx=1:nvar
			if (isnumeric(afun))
				a = afun(idx,jdx);
			else
				a  = afun(idx,jdx,n);
			end
			di(:,:,4+jdx,idx) = afun(idx,jdx,n);
		end
		di(:,:,4+idx,idx) = di(:,:,4+idx,idx) + dt*2*(d(1,idx)/dx(1)^2+d(2,idx)/dx(2)^2);
	end
	%di = reshape(di,n(1),n(2),4+nvar,nvar);
end % function


