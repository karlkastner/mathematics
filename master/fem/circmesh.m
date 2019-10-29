% Wed Feb 27 00:32:57 MSK 2013
% Karl Kästner, Berlin

% generates a circular mesh
function [P T B X] = circmesh(n, r)
	[P T B X] = mesh_2d_uniform([31 31], [r r], 0.5*[r r]);
	%X = ones(n,1)*(0:n-1) * 2/(n-1) -1;
	%Y = X';
	%X = X(:);
	%Y = Y(:);

%	Rs = sum(P.^2,2);
%	Ip = (Rs) < 0.25*(r + r/n)^2;
%	It = Ip(T(:,1)).*Ip(T(:,2)).*Ip(T(:,3));
%	tdx = find(It);
%	T = T(tdx,:);

	X = P(:,1);
	Y = P(:,2);
	X = X+eps;
	Y = Y+eps;
	R = r./min(sqrt(1+(X./Y).^2), sqrt(1+(Y./X).^2));
	X = X.*R;
	Y = Y.*R;
	P = [X Y];
end

