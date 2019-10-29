% Thu Apr 19 13:15:39 MSK 2012
% Karl Kästner, Berlin

% generate mesh
function [P T Bc] = mesh_1d_uniform(n, L)
	% points
	if (nargin == 1)
		P = n;
		n = length(n);
	else
		P = L*(0:n-1)'/(n-1);
	end
	% line elements
	T = [ (1:n-1)' (2:n)' ];
	% boundary points
	Bc = [1 1;
              n 2];
end

