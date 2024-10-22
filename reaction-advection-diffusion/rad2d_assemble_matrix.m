% 2024-01-31 11:57:44.706345420 +0100
% Karl Kastner, Berlin
%
%% assemble the discretization matrix of a two-dimensional 
%% reaction-advection-diffusion equation
function A = rad2d_assemble_matrix(a,ad,e,n,L,nvar)
	dx = L./n;
	%I  = speye(prod(n));
	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,dx,2,{'circular','circular'},true);
	%D2 = D2x+D2y;
	A = [];
	for idx=1:nvar
		Ai = [];
	for jdx=1:nvar
		if (idx==jdx)
			if (iscell(a))
				Ai = [Ai,diag(sparse(flat(a{idx,jdx}))) ...
				         + e(1,idx)*D2x ...
				         + e(2,idx)*D2y ...
					 ... % + kernel2matrix(ad(:,:,idx),n,dx) ...
					];
			else
				Ai = [Ai,diag(sparse(squeeze(a(idx,jdx,:)))) ...
					+ e(1,idx)*D2x ... 
					+ e(2,idx)*D2y ... 
					... % + kernel2matrix(ad(:,:,idx),n,dx) ...
					];
			end
		else
			if (iscell(a))
				Ai = [Ai,diag(sparse(flat(a{idx,jdx})))];
			else
				Ai = [Ai,diag(sparse(squeeze(a(idx,jdx,:))))];
			end
		end % if
	%A  = [diag(sparse(flat(a{1,1}))) - d(1)*D2, diag(sparse(flat(a{1,2})));
        %      diag(sparse(flat(a{2,1}))), diag(sparse(flat(a{2,2}))) - d(2)*D2];
	end % jdx
		A  = [A;Ai];
	end % idx
end % assemble_rad_nvar

