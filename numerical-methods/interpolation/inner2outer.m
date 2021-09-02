% Fri  1 Jun 10:34:53 CEST 2018
%% linear interpolation of segment mit point to grid points at segment ends
%% assumes equal grid spacing
function x = inner2outer(x,dim)
	if (isvector(x))
		x=cvec(x);
	end
	if (nargin()<2)
		dim = 1;
	end
	if (size(x,1) == 1)
		x = [x;x];
	else
		if (1 ~= dim)
			x = shiftdim(x,dim-1);
		end
	        x = [1.5*x(1,:,:,:) - 0.5*x(2,:,:,:);
                     mid(x);
	             1.5*x(end,:,:,:) - 0.5*x(end-1,:,:,:);
                    ];
		if (1 ~= dim)
			nd = ndims(x);
			x = shiftdim(x,nd-dim+1);
		end
	end
end

