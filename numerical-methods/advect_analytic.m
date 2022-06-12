% Mon 30 May 13:23:07 CEST 2022
% analytic solution to the advection equation
% solve advection equation in 1 or two dimensions
% with constant coefficient and circular boundary conditions
function y = advect_analytic(t,y0,L,a)
	if (isvector(y0))
		y0 = cvec(y0);
	end
	y = advect(t,y0,L(1),a(1));
	% second dimension
	if (size(y,2) > 1)
		y = advect(t,y.',L(2),a(2)).';
	end
end

function y = advect(t,y0,L,a)
	nx = size(y0,1);
	dx = L./nx;
	% this is just a linear interpolation
	dxt  = t*a/dx;
	dxt_ = floor(dxt);
	q    = (dxt - dxt_);

	idl = circshift(1:nx,-dxt_);
	idr = circshift(1:nx,-dxt_-1);

	y = (1-q)*y0(idl,:) + q*y0(idr,:);
end

