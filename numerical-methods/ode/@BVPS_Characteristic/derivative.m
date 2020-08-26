	% derivative of solution
	if (nargout > 3)
	% TODO, this is not working any more, due to regression
	% since first and second order odes were coupled
	l = [];
	m = 0;
	dydx    = [ (  l(:,1).*ypm(1:m:m*nxc-m+1).*obj.exp(-l(:,1).*dx(1:nxc)/2) ...
                     + l(:,2).*ypm(3:m:m*nxc    ).*obj.exp(-l(:,2).*dx(1:nxc)/2));
	            (  l(nxc,1).*ypm(m*nxc-2)*obj.exp(l(nxc,1)*dx(nxc)/2) ...
                     + l(nxc,2).*ypm(  m*nxc)*obj.exp(l(nxc,2)*dx(nxc)/2)) ];
	end

