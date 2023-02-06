% Fri 15 Dec 11:00:30 CET 2017
% Karl Kastner, Berlin
%
%% first derivative on variable mesh
%% second order accurate
% function dy_dx = derivative1(x,y,order)
function dy_dx = derivative1(x,y,order)
	istransposed = false;
	if (isvector(y) && isrow(y))
		y = cvec(y);
		istransposed = true;
	end

	D1 = derivative_matrix_1_1d(x,[]);

	dy_dx = D1*y;

	if (istransposed)
		dy_dx = dy_dx.';
	end
end

