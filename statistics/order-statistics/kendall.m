% Sa 1. Aug 12:48:22 CEST 2015
% Karl Kastner, Berlin
%
%% kendall correlation coefficient
%
function rho = kendall(x,y)
	if (isvector(x))
		x = cvec(x);
		y = cvec(y);
	end

	n  = size(x,1);
	m  = size(x,2);
	S  = zeros(1,m);
	N  = zeros(1,m);
	tx = zeros(1,m);
	ty = zeros(1,m);
	for idx=1:n
	 for jdx=idx+1:n
		sx = sign(x(idx,:)-x(jdx,:));
		sy = sign(y(idx,:)-y(jdx,:));
		S  = S+sx.*sy;
		tx = tx + (1-abs(sx)).*abs(sy);
		ty = ty + (1-abs(sy)).*abs(sx);
	        N  = N+abs(sx).*abs(sy);
	 end % for jdx
	end % for idx
	% number of ties
	rho = S./sqrt((N+tx).*(N+ty));
%	rho = S*1/(n*(n-1));
end % kendall

