% 2024-01-24 13:09:21.118254769 +0100
function jacobi_step(obj,k)

	obj.resfun(k);

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

end % jacobi_step

