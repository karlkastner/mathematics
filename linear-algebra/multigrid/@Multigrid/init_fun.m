% Tue 24 Sep 18:38:49 CEST 2024
% Karl Kastner, Berlin
%
%% init multigrid, determine matrices at coarser levels through a function
function init_fun(obj,fun,L,n,nvar,xi)
	obj.L    = L;
	obj.n    = n;
	obj.nvar = nvar;
	obj.fun  = fun;
	%xi = reshape(xi,n(1)*n(2),nvar);
	obj.s(1).diagonals = fun(L,n,nvar,xi);
	obj.s(1).dx = L./n(1:2);
	obj.s(1).res = zeros(n(1),n(2),obj.nvar);
	k = 1;
	while (n(1)>1)
		k = k+1;
		n = n/2;
		if (~isempty(xi))
		% note that if xi can have more slices than nvar in case of optional arguments
		n3  = size(xi,3);
		xi_ = zeros(n(1),n(2),n3);
		for idx=1:n3
			xi_(:,:,idx) = downsample_2d(xi(:,:,idx),'fw');
		end
		end
		xi = xi_;
		% diagonals [n(1)*n(2), 4+nvar, nvar]
		obj.s(k).diagonals = fun(L,n,nvar,xi);
		obj.s(k).res = zeros(n(1),n(2),obj.nvar);
		obj.s(k).x   = zeros(n(1),n(2),obj.nvar);
	end % while
end % init

