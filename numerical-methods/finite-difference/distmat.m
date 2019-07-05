% 2013-08-07 00:47:18.000000000 +0200
% Karl Kastner, Berlin
%
%% distance matrix for a 2 dimensional rectangular matrix
%
% TODO use derivative matrices
function W = distmat(Z,flag)
	ny = size(Z,1);
	nx = size(Z,2);
	s = size(Z);
	m = 0;
	switch (flag)
	% 4-mneighbourhood
	case {4,'4'}
	for idx=1:ny
	 for jdx=1:nx
		if (idx-1 > 0)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx-1,jdx), (Z(idx,jdx) - Z(idx-1,jdx)).^2];
		end
		if (idx + 1 <= ny)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx+1,jdx), (Z(idx,jdx) - Z(idx+1,jdx)).^2];
		end
		if (jdx-1 > 0)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx,jdx-1), (Z(idx,jdx) - Z(idx,jdx-1)).^2];
		end
		if (jdx+1 <= nx)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx,jdx+1), (Z(idx,jdx) - Z(idx,jdx+1)).^2];
		end
	 end % jdx
	end % idx
	% 8 neighbourhood
	case {8, '8'}
	for idx=1:ny
	 for jdx=1:nx
		if (idx-1 > 0 && jdx-1>0)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx-1,jdx-1), (Z(idx,jdx) - Z(idx-1,jdx-1)).^2];
		end
		if (idx-1 > 0)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx-1,  jdx), (Z(idx,jdx) - Z(idx-1,jdx)).^2];
		end
		if (idx-1 > 0 && jdx+1<=nx)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx-1,jdx+1), (Z(idx,jdx) - Z(idx-1,jdx+1)).^2];
		end

		if (jdx-1 > 0)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,  idx,jdx-1), (Z(idx,jdx) - Z(idx,jdx-1)).^2];
		end
		if (jdx+1 <= ny)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx,  jdx+1), (Z(idx,jdx) - Z(idx,jdx+1)).^2];
		end
		if (idx+1 <= ny && jdx-1 > 0)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx+1,jdx-1), (Z(idx,jdx) - Z(idx+1,jdx-1)).^2];
		end
		if (idx+1 <= ny)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx+1,  jdx), (Z(idx,jdx) - Z(idx+1,jdx)).^2];
		end
		if (idx+1 <= ny && jdx+1 <= nx)
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx), sub2ind(s,idx+1,jdx+1), (Z(idx,jdx) - Z(idx+1,jdx+1)).^2];
		end
	 end % jdx
	end % idx
	case {inf,'inf'}
		for idx=1:ny
		 for jdx=1:nx
		  for kdx=1:ny
                   for ldx=1:nx
			m = m+1;
			buf(m,1:3) = [sub2ind(s,idx,jdx) sub2ind(s,kdx,ldx) 0.5*(Z(idx,jdx) - Z(kdx,ldx)).^2 ];
		   end % ldx
		  end % kdx
		 end % jdx
		end % idx
	case {'rand'}
		for idx=1:ny
		 for jdx=1:nx
		  for mdx=1:10
			kdx=randi(ny);
			ldx=randi(nx);
			m = m+1;
			% TODO this causes problems, when indices are chosen twice, as sparse adds duplicates
			buf(m,1:3) = [sub2ind(s,idx,jdx) sub2ind(s,kdx,ldx) 0.5*(Z(idx,jdx) - Z(kdx,ldx)).^2 ];
			end
			end
		end	
	end % switch
	W = sparse(buf(:,1),buf(:,2),buf(:,3),ny*nx,ny*nx);
	% TODO, that is a quick hack to make the mat symmetic, this doubles certain distances
	W = W + W';
%	this does not work, as transposing is not equal to swapping x and y
%	ind1  = find(W(:) > 0);
%	[y x] = ind2sub(size(W),ind1);
%	ind2  = sub2ind(size(W),x,y);
%	W(ind2) = W(ind1);
%	W(y,x) = W(x,y);
	

%	condest(W)
%	buf(fdx,3) = 1./buf(fdx,3);
%	W = W+W';	% todo, exploit symmetry
end % function distmat

