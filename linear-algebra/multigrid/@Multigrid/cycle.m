% 2024-01-24 12:40:44.370170282 +0100
function cycle(obj,k)
	if (size(obj.s(k).x,1) > 1)
		% pre-smooth, hackbusch 2.5.2.b
		for idx=1:obj.m
			obj.jacobi_step(k);
		end	

		for gdx=1:obj.gamma
		% n.b. for the w-cycle, this is simply repated twice
		% residual
		obj.resfun(k);

		% downsample (restrict)
		for idx=1:obj.nvar
			obj.s(k+1).b(:,:,idx) = downsample_2d(obj.s(k).res(:,:,idx));
		end

		% initialize
		obj.s(k+1).x(:)=0;
	
		% recurse to approximate error
		obj.cycle(k+1);

		% upsample and correct
		for idx=1:obj.nvar
			obj.s(k).x(:,:,idx) = obj.s(k).x(:,:,idx) - upsample_2d(obj.s(k+1).x(:,:,idx));
		end

		end % for gamma

		% post-smooth
		for idx=1:obj.m
			obj.jacobi_step(k);
		end
	else
		% solve the nvxnv system here
		A = cell2mat(obj.s(k).a);
if (1)
		obj.s(k).x(1,1,1:obj.nvar) = A \ cvec(squeeze(obj.s(k).b));
else
		obj.s(k).x(1,1,1:obj.nvar) = cvec(squeeze(obj.s(k).b)) ./ diag(A); 
end
%		for idx=1:obj.nvar
%			obj.s(k).x(idx) = obj.s(k).b(idx)./(obj.s(k).a(:,:,idx));
%			obj.s(k).x(idx) = obj.s(k).b(idx)./(obj.s(k).a(:,:,idx));
%		end	
	end
end % cycle	

