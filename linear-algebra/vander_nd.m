% Sun 13 Feb 11:47:51 CET 2022
% n-dimensional vandermonde matrix
% x : columns : variables (x,y,z,...), rows: samples (x1;x2;x3 ...)
function [A,id] = vander_nd(x,order)
	n = size(x);
	switch (order)
	case {0}
		A = ones(n(1),1);
		id = [0;0];
	case {1}
		A = [ones(n(1),1),x];
		id = [[0;0],[(1:n(2));zeros(n(2),1)]];
	case {2}
		% constant and linear term
		A  = [ones(n(1),1),x,zeros(n(1),0.5*n(2)*(n(2)+1))];
		id = [[0;0],[(1:n(2));zeros(1,n(2))],zeros(2,0.5*n(2)*(n(2)+1))];
		k  = n(2)+1;
		% quadratic terms
		for idx=1:n(2)
			for jdx=idx:n(2)
				k = k+1;
				id(1,k) = idx;
				id(2,k) = jdx; 
				A(:,k)  = x(:,idx).*x(:,jdx);
			end
 		end
	otherwise
		error('not yet implemented');
	end
end

