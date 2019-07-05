% Fri 15 Dec 11:00:30 CET 2017
% Karl Kastner, Berlin
% 
%% second derivative on a variable mesh
function d2y_dx2 = derivative2(x,y,order)

	transposed = false;
	if (isvector(y) && isrow(y))
		y = cvec(y);
		transposed = true;
	end
	%dy = zeros(size(y));

	D2  = derivative_matrix_2_1d(x);
	d2y_dx2 = D2*y;

	if (transposed)
		% retranspose
		d2y_dx2 = d2y_dx2.';
	end
	%,L,order,bcl,bcr)

%	for idx=1:size(y,2)
%		dy(1,idx) = (y(2,idx)-y(1,idx))/(x(2)-x(1));
%		dy(2:end-1,idx) = [	        (x(2:end-1) - x(3:end))./((x(1:end-2) - x(2:end-1)).*(x(1:end-2) - x(3:end))).*y(1:end-2,idx) ...
%				  + (x(1:end-2) - 2*x(2:end-1) + x(3:end))./((x(1:end-2) - x(2:end-1)).*(x(2:end-1) - x(3:end))).*y(2:end-1,idx) ...
%				  +      (x(2:end-1) - x(1:end-2))./((x(1:end-2) - x(3:end)).*(x(2:end-1) - x(3:end))).*y(3:end,idx) ];
%		dy(end,idx) = (y(end)-y(end-1))./(x(end)-x(end-1));
%	end % for idx
end

