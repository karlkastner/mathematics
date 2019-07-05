% Sa 6. Feb 23:40:22 CET 2016
%
%% eigenvalues of stacked 2x2 matrices
function [V, E] = eig2x2(a,b,c,d)
	if (1 == nargin())
		d = shiftdim(a(2,2,:),+1);
		c = shiftdim(a(2,1,:),+1);
		b = shiftdim(a(1,2,:),+1);
		% must come last as a is overwritten
		a = shiftdim(a(1,1,:),+1);
	end
	n = length(a);

	% determinant
	det = a.*d - b.*c;
	% trace
	trace = a+d;
	% eigenvalues
	e1 = 0.5*trace + sqrt(0.25*trace.^2-det);
	e2 = 0.5*trace - sqrt(0.25*trace.^2-det);
	% stack here, to avoid empty matrix 0x1 vs 0 scalar missmatch
	E = [e1; e2];
	if (nargout() > 1)
		V   = ones(2,2,n);

		% TODO not exactly zero, but when term under sqrt almost zero
		fdx        = (c==0);
		V(2,1,fdx) = 0;
		V(1,2,fdx) = b(fdx)./(d(fdx)-a(fdx));
		%V(1,2,fdx) = -d(fdx)./(a(fdx)-d(fdx));
		fdx        = ~fdx;
		V(1,1,fdx) = (e1(fdx)-d(fdx))./c(fdx);
		V(1,2,fdx) = (e2(fdx)-d(fdx))./c(fdx);

		% normalize
		n1       = sqrt(V(1,1,:).^2 + V(2,1,:).^2);
		V(1,1,:) = V(1,1,:)./n1;
		V(2,1,:) = V(2,1,:)./n1;

		n2       = sqrt(V(1,2,:).^2 + V(2,2,:).^2);
		V(1,2,:) = V(1,2,:)./n2;
		V(2,2,:) = V(2,2,:)./n2;

%	% eigenvector norms
%	% this is for avoidance of division by zero, but also should improve conditioning
%	% TODO, best choice for condition?
%	fdx        = abs(c) + abs(d) > abs(a) + abs(b);
%	
%	n1         = 1./sqrt(c(fdx).^2 + (e1(fdx)-d(fdx)).^2)
%	V(1,1,fdx) = (e1(fdx)-d(fdx)).*n1;
%	V(2,1,fdx) = c(fdx).*n1;
%	n2         = 1./sqrt(c(fdx).^2 + (e2(fdx)-d(fdx)).^2)
%	V(1,2,fdx) = (e2(fdx)-d(fdx)).*n2;
%	V(2,2,fdx) = c(fdx).*n2;
%
%	fdx        = ~fdx;
%	n1         = 1./sqrt(b(fdx).^2 + (e1(fdx)-a(fdx)).^2)
%	V(1,1,fdx) = b(fdx).*n1;
%	V(2,1,fdx) = (e1(fdx)-a(fdx)).*n1;
%	n1         = 1./sqrt(b(fdx).^2 + (e2(fdx)-a(fdx)).^2)
%	V(1,2,fdx) = b(fdx).*n1;
%	V(2,2,fdx) = (e2(fdx)-a(fdx)).*n1;
	else
	V = E;
	end
end

