% Sat Feb 25 17:54:38 MSK 2012
% Karl Kästner, Berlin

function [P T B] = mesh_3d_uniform(n, L, x0)
	if (nargin() < 3)
		x0 = [0 0 0];
	end

	% minimal meshes
	if (-1 == n(1))
		P = [0 0 0;
                     1 0 0;
                     0 1 0;
                     0 0 1];
		T = [1 2 3 4];
		B = [1 2 3;
                     1 2 4;
                     1 3 4;
                     2 3 4];
		return;
	end
	if (-2 == n(1))
		P = [0 0 0;
                     1 0 0;
                     0 1 0;
                     0 0 1;
                    -1 0 0];
		P(:,1) = 0.5*P(:,1)+0.5;
		T = [2 4 1 3;
		     1 5 3 4];
		B = [1 2 3;
                     1 2 4;
                     2 3 4;
		     1 5 3;
                     1 5 4;
                     3 5 4];
		return;
	end

	% point coordinates per axis
	X1 = L(1)*(0:n(1)-1)'/(n(1)-1);
	X2 = L(2)*(0:n(2)-1)'/(n(2)-1);
	X3 = L(3)*(0:n(3)-1)'/(n(3)-1);

	% shift
	X1 = X1 - x0(1);
	X2 = X2 - x0(2);
	X3 = X3 - x0(3);

	% 3D mesh
	I1 = ones(n(1),1);
	I2 = ones(n(2),1);
	I3 = ones(n(3),1);
	P = [ kron(kron(X1,I2),I3), kron(kron(I1,X2),I3), kron(kron(I1,I2),X3) ];

	% tetrahedra
	tdx=0;
	T = zeros(6*prod(n-1), 4);

	% TODO generate boundaries correctly for different number of points per axis
	%n = n(1);

	% boundary faces
	if (nargout() > 2)
	bdx=0;
	B = zeros(4*((n(1)-1)*(n(2)-1) + (n(1)-1)*(n(3)-1) + (n(2)-1)*(n(3)-1)), 3);
	for xdx=1:n(1)-1
	  for ydx=1:n(2)-1
	    for zdx=1:n(3)-1
		% point index
		k = xdx + (ydx-1)*n(2) + (zdx-1)*n(3)^2;
		% local triangluation of a rectangular box into six tetraeders
		T(tdx+1:tdx+6,:) = [   k+0,       k+1,     k+n(1)+1,   k+n(1)^2;
				       k+1,   k+n(1)^2+1, k+n(1)^2+n(1)+1,   k+n(1)^2;
				       k+1,     k+n(1)+1, k+n(1)^2+n(1)+1,   k+n(1)^2;
				       k+0,     k+n(1)+1,       k+n(1), k+n(1)^2+n(1);
				       k+0,     k+n(1)+1,   k+n(1)^2+n(1),   k+n(1)^2;
				     k+n(1)+1, k+n(1)^2+n(1)+1,   k+n(1)^2+n(1),   k+n(1)^2];
		tdx = tdx+6;
		if (1==zdx)
			% top faces
			bdx=bdx+1; B(bdx,:) = [k+0 k+1   k+n(1)+1];
			bdx=bdx+1; B(bdx,:) = [k+0 k+n(1)+1 k+n(1)];
		end
		if (n(3)-1==zdx)
			% bottom faces
			bdx=bdx+1; B(bdx,:) = [k+n(1)^2+1   k+n(1)^2+n(1)+1 k+n(1)^2];
			bdx=bdx+1; B(bdx,:) = [k+n(1)^2+n(1)+1 k+n(1)^2+n(1)   k+n(1)^2];
		end
		if (1==xdx)
			% left faces
			bdx=bdx+1; B(bdx,:) = [k+0 k+n(1)^2 k+n(1)^2+n(1)];
			bdx=bdx+1; B(bdx,:) = [k+0 k+n(1)   k+n(1)^2+n(1)];
		end
		if (n(1)-1==xdx)
			% right faces
			bdx=bdx+1; B(bdx,:) = [k+1 k+n(1)^2+1 k+n(1)^2+n(1)+1];
			bdx=bdx+1; B(bdx,:) = [k+1 k+n(1)+1   k+n(1)^2+n(1)+1];
		end
		if (1==ydx)
			% back faces
			bdx=bdx+1; B(bdx,:) = [k+0 k+1     k+n(1)^2];
			bdx=bdx+1; B(bdx,:) = [k+1 k+n(1)^2 k+n(1)^2+1];
		end
		if (n(2)-1==ydx)
			% fron(1)t faces
			bdx=bdx+1; B(bdx,:) = [k+n(1)   k+n(1)+1     k+n(1)^2+n(1)];
			bdx=bdx+1; B(bdx,:) = [k+n(1)+1 k+n(1)^2+n(1) k+n(1)^2+n(1)+1];
		end
	    end % zdx
	  end % ydx
	end % xdx
	end
end % triangulate_3d

