% 2024-01-24 13:09:21.118254769 +0100
%% single jacobi iteration
function jacobi_step(obj,k)

	obj.resfun(k);

	% TODO this misses the remaining reaction terms
	if (isempty(obj.fun))
	   for idx=1:obj.nvar
			% TODO, why no - in front of a?
			d = (    obj.s(k).a{idx,idx} ...
			    -    obj.ad(1,2,idx)/obj.s(k).dx(1) ...
			    + 2*obj.e(1,idx)/obj.s(k).dx(1)^2 ... 
			    -    obj.ad(2,2,idx)/obj.s(k).dx(2) ...
			    + 2*obj.e(2,idx)/obj.s(k).dx(2)^2  ) ./ obj.o;
			obj.s(k).x(:,:,idx) = ( obj.s(k).x(:,:,idx) ...
				- (obj.s(k).res(:,:,idx)./d) );
	   end % for idx
	
	else % of if isempty obj.fun

	% A*x = rhs
	% (D+R)*x = rhs
	% x = (rhs - R*x)*D
	% x = (rhs - R*x - D*x + D*x)/D
	% x = (rhs - A*x + D*x)/D
	% x = x + (rhs - A*x)/D
	% x = x + (-res)/D
	    for idx=1:obj.nvar
		obj.s(k).x(:,:,idx) = (    obj.s(k).x(:,:,idx) ...
 			                - ( obj.s(k).res(:,:,idx) ...
					   ./ (1./obj.o*obj.s(k).diagonals(:,:,5)) ...
					  ) ...
				      );
	    end % for idx
	end % of else of if isempty obj.fun
end % jacobi_step

