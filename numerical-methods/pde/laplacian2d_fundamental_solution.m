% Wed  5 Feb 20:12:50 +08 2020
function phi = laplacian2d_fundamental_solution(x0,y0)
	if (nargin()<1)
		x0 = 0;
		y0 = 0;
	end
	syms x y
	phi = -1/sym(pi)*log(sqrt((x-x0).^2 + (y-y0).^2));
end
