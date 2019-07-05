% Thu 11 Jan 10:15:43 CET 2018
% Karl Kastner, Berlin
%
%% flag points contained in tetrahedron
%
function [in] = inTetra(X,Y,Z,X0,Y0,Z0)
	n = size(X,1);
	c = zeros(n,1);

	% translate so that P1 is at origin
	% to reduce round off error
	P2  = [X(:,2) - X(:,1), Y(:,2) - Y(:,1), Z(:,2) - Z(:,1)].';
	P3  = [X(:,3) - X(:,1), Y(:,3) - Y(:,1), Z(:,3) - Z(:,1)].';
	P4  = [X(:,4) - X(:,1), Y(:,4) - Y(:,1), Z(:,4) - Z(:,1)].';
	P0  = [X0(:)  - X(:,1), Y0(:)  - Y(:,1), Z0(:)  - Z(:,1)].';
	P1  = zeros(size(P2));


	% a point is inside an n-simplex if it is on the
	% same side of all (n-1)-simplexes,
	% i.e. behind all walls in 3d, behind all lines in 2d,
	% on the inside side of the vertices in 1d

	v0 = vander_3d(P1,P2,P3,P4,1);
	d0 = det4x4(v0);
	v1 = vander_3d(P0,P2,P3,P4,1);
	d1 = det4x4(v1);
	v2 = vander_3d(P1,P0,P3,P4,1);
	d2 = det4x4(v2);
	v3 = vander_3d(P1,P2,P0,P4,1);
	d3 = det4x4(v3);
	v4 = vander_3d(P1,P2,P3,P0,1);
	d4 = det4x4(v4);
	
	s0 = sign(d0);
	s1 = sign(d1);
	s2 = sign(d2);
	s3 = sign(d3);
	s4 = sign(d4);

	in =   (s0==sign(d1) | 0==sign(d1)) ...
	     & (s0==sign(d2) | 0==sign(d2)) ...
	     & (s0==sign(d3) | 0==sign(d3)) ...
	     & (s0==sign(d4) | 0==sign(d4));
	
	in = squeeze(in);
	in = cvec(in);

end % inTetra

