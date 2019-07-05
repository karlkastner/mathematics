% Mo 2. Nov 18:40:31 CET 2015
% Karl Kastner, Berlin
% centre and radius of circle from 3 points
%% excircle of a triangle
function [Xc Yc R] = tri_excircle(X,Y)
	n = size(X,1);
	R  = [];

	if (~issym(X))	
		A = zeros(2,2,n);
		b = zeros(2,n);
	end

	A(1,1,:) = -(Y(:,1)-Y(:,2));
	A(1,2,:) = -(Y(:,3)-Y(:,2));
	A(2,1,:) =   X(:,1)-X(:,2);
	A(2,2,:) =   X(:,3)-X(:,2);
	b(1,:)   = 0.5*(X(:,1)-X(:,3));
	b(2,:)   = 0.5*(Y(:,1)-Y(:,3));

	% solve equation systems
	Ai = inv2x2(A);
	%s  = matvec2x2(Ai,b);
	s  = matvec3(Ai,b);

	Xc(:,1) = 0.5*(X(:,1)+X(:,2)) + s(1,:).'.*(Y(:,1)-Y(:,2));
	Yc(:,1) = 0.5*(Y(:,1)+Y(:,2)) - s(1,:).'.*(X(:,1)-X(:,2));

	if (n > 0 && nargout > 2)
		if (~issym(X))
			R = hypot(X(:,1)-Xc, Y(:,1)-Yc);
		else
			R = sqrt( (X(:,1)-Xc).^2 + (Y(:,1)-Yc).^2 );
		end
	end
end % circumferencecircle

