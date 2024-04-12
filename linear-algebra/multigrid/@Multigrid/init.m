% 2024-01-10 19:59:21.792153268 +0100
function init(obj,a,ad,e,L,n,nvar)
	if (nargin()>5)
		nvar = 1;
	end
	if (isempty(ad))
		ad = zeros(2,3,nvar);
	end
	obj.L = L;
	obj.n = n;
	obj.ad = ad;
	obj.e = e;
	obj.nvar = nvar;
	obj.s(1).a = a;
	obj.s(1).dx = L./n;
	obj.s(1).res = zeros(n(1),n(2),obj.nvar);
	k = 1;
	while (n(1)>1)
		k=k+1;
		n=n/2;
		for idx=1:obj.nvar
		    for jdx=1:obj.nvar
			if (isscalar(a{idx,jdx}))
				obj.s(k).a{idx,jdx} = a{idx,jdx};
			else
				obj.s(k).a{idx,jdx} = downsample_2d(obj.s(k-1).a{idx,jdx});
			end
			%obj.s(k).a(:,:,nvar) = downsample_2d(obj.s(k-1).a(:,:,nvar));
		    end % for jdx
		end % for idx
		obj.s(k).dx  = L./n;
		obj.s(k).res = zeros(n(1),n(2),obj.nvar);
		obj.s(k).x   = zeros(n(1),n(2),obj.nvar);
	end % while
end % init

