% Thu May 17 14:23:03 MSK 2012
% Karl Kästner
%
%function V = vander_2d(X,Y,p)
%
% set up the vandermonde matrix in two dimensions
%
function V = vander_2d(X,Y,p)
%	if (nargin() < 3)
%		n = py;
%		py = px(:,2);
%		px = px(:,1);
%	end
	siz = size(X);
	Xp = ones(siz(1),siz(2),p+1);
	Yp = ones(siz(1),siz(2),p+1);
	if (issym(X))
		Xp = sym(Xp);
		Yp = sym(Yp);
	end

	% powers of x and y
	for jdx=2:p+1
		Xp(:,:,jdx) = X.*Xp(:,:,jdx-1);
		Yp(:,:,jdx) = Y.*Yp(:,:,jdx-1);
	end
	

	% preallocate memory
	n = (p+1)*(p+2)/2;
	V = zeros(size(Xp,1),n,size(Xp,2));
	if (issym(X))
		V = sym(V);
	end
	k=0;
	for pdx=0:p
		for xp=pdx:-1:0
			% power of y
			yp = pdx-xp;
			k = k+1;
			V(:,k,:) = Xp(:,:,xp+1).*Yp(:,:,yp+1);
		end % for xp
	end % for pdx
end % function vander_2d

