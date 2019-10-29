% Thu May 17 14:23:03 MSK 2012
% Karl Kästner

% set up the vandermonde matrix in two dimensions
function V = vander_2d(p,n)
	% preallocate memory
	if (isnumeric(p))
		V = zeros(size(p,1),0.5*(n+1)*(n+2));
	end
	m=1;
	for idx=0:n
		for jdx=0:idx
			% todo, improve without transcendental functions
			V(:,m) = p(:,1).^(idx-jdx).*p(:,2).^(jdx);
			m = m+1;
		end % for jdx
	end % for idx
end % function vander_2d

