% 2024-01-10 19:59:21.792153268 +0100
% Karl Kastner, Berlin
%
%% set up multigrid, and matrix at coarser levels
%
%% number of coefficients    : (4 + nvar^2)*n(1)*n(2)
%% number of diagonals       : ndiag = (4 + nvar)
%% number of matrix elements : nnz   = n(1)*n(2)*(4+nvar)*nvar
function init(obj,a,ad,e,L,n,nvar)
	if (nargin()>5)
		nvar = 1;
	end
	if (isempty(ad))
		ad = zeros(2,3,nvar);
	end
	obj.L    = L;
	obj.n    = n;
	obj.ad   = ad;
	obj.e    = e;
	obj.nvar = nvar;
	obj.s(1).a = a;
	obj.s(1).dx = L./n;
	obj.s(1).res = zeros(n(1),n(2),obj.nvar);
	k        = 1;
	while (n(1)>1)
		k=k+1;
		n=n/2;
%		if (obj.opt.operator_dependent_grid_transfer)
%		for idx=1:obj.nvar
%			% downsample the diffusion coefficients
%			obj.s(k).ex{idx} = downsample_2d(obj.s(k-1).ex,obj.opt.downsampling_mode);
%			obj.s(k).ey{idx} = downsample_2d(obj.s(k-1).ey,obj.opt.downsampling_mode);
%			% downsample the advection coefficients
%			obj.s(k).ax{idx} = downsample_2d(obj.s(k-1).ax,obj.opt.downsampling_mode);
%			obj.s(k).ay{idx} = downsample_2d(obj.s(k-1).ay,obj.opt.downsampling_mode);
%		end % for idx
%		end % if odgt
		% downsample the recation coefficients
		for idx=1:obj.nvar
		    for jdx=1:obj.nvar
			if (isscalar(a{idx,jdx}))
				obj.s(k).a{idx,jdx} = a{idx,jdx};
			else
				obj.s(k).a{idx,jdx} = downsample_2d(obj.s(k-1).a{idx,jdx},obj.opt.downsampling_mode);
			end
			%obj.s(k).a(:,:,nvar) = downsample_2d(obj.s(k-1).a(:,:,nvar));
		    end % for jdx
		end % for idx
		obj.s(k).dx  = L./n;
		obj.s(k).res = zeros(n(1),n(2),obj.nvar);
		obj.s(k).x   = zeros(n(1),n(2),obj.nvar);
		% set up the coefficients
%		if (obj.opt.operator_dependent_grid_transfer)
%			obj.s(k).diags = rad2d_coeffs2grid(obj.s(k).ex,obj.s(k).ey,obj.s(k).ax,obj.s(k).ay,obj.s(k).a,obj.s(k).dx);
%		end % if odgt
	end % while
end % init

