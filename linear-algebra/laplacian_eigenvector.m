% 2022-06-13 10:05:51.516852870 +0200
% TODO, different BCS
% 2d vector is v = kron(v1,v2) 
function [v1,v2] = laplacian_eigenvector(n,id)
	x = (1:n(1))';
	v1 = sin(pi*id(1)/(n(1)+1)*x);
	v1 = sqrt(2/(n(1)+1))*v1;
	if (length(n) == 2)
		x = (1:n(2))';
		v2 = sin(pi*id(2)/(n(2)+1)*x);
		v2 = sqrt(2/(n(2)+1))*v2;
	end
end

