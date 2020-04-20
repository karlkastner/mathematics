% Fri 15 Dec 11:00:30 CET 2017
% Karl Kastner, Berlin
%
%% first derivative on variable mesh
%% second order accurate
function dy_dx = derivative1(x,y,order)
	istransposed = false;
	if (isvector(y) && isrow(y))
		y = cvec(y);
		istransposed = true;
	end

	D1 = derivative_matrix_1_1d(x,[]);

	dy_dx = D1*y;

%	dy_dx = zeros(size(y));
%	for idx=1:size(y,2)
%		dy_dx(1,idx) = (y(2,idx)-y(1,idx))/(x(2)-x(1));
%		dy_dx(2:end-1,idx) = [	        (x(2:end-1) - x(3:end))./((x(1:end-2) - x(2:end-1)).*(x(1:end-2) - x(3:end))).*y(1:end-2,idx) ...
%				  + (x(1:end-2) - 2*x(2:end-1) + x(3:end))./((x(1:end-2) - x(2:end-1)).*(x(2:end-1) - x(3:end))).*y(2:end-1,idx) ...
%				  +      (x(2:end-1) - x(1:end-2))./((x(1:end-2) - x(3:end)).*(x(2:end-1) - x(3:end))).*y(3:end,idx) ];
%		dy_dx(end,idx) = (y(end)-y(end-1))./(x(end)-x(end-1));
%	end % for idx

	if (istransposed)
		dy_dx = dy_dx.';
	end
end

