function connectivity(A)
	s = size(A);
	buf = [];
	for idx=1:s(1)
	for jdx=1:s(2)
		% for each neighbour
		for dx=-1:1
		for dy=-1:1
			% TODO boundary
			dz = A(idx,jdx) - A(idx+dx,jdx+dy);
			k = k+1;
			buf(k,2) = idx+(jdx-1)*s(1);
			buf(k,2) = (idx+dx) + (jdx+dy-1)*s(2);
			buf(k,3) = sqrt(dx^2 + dy^2 + dz^2);
		end
		end
	end
	end
end
