% Fri 21 Aug 12:02:13 +08 2020
function dy_u  = upwind_difference(x,y)
	dy   = diff(y);
	% heed x-direction
	s    = sign(x(2)-x(1));
	% flag for upwind direction
	fdx  = s.*y(2:end-1) < 0;
	dy_u = fdx.*dy(2:end) + (1-fdx).*dy(1:end-1);
end

