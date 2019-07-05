% 2015-11-02 18:56:47.018003787 +0100
%% cartesian to barycentric coordinates
function c = tobarycentric2(P1,P2,P3,P0)
	n = size(P0,1);
	c = zeros(3,n);

	% translate so that P1 is at origin
	P2 = P2-P1;
	P3 = P3-P1;
	P0 = P0-P1;

	if (1)
	A = zeros(2,2,n);
	A(1,1,:) = P2(:,1);
	A(1,2,:) = P3(:,1);
	A(2,1,:) = P2(:,2);
	A(2,2,:) = P3(:,2);
	b        = P0';
	else

	% because X1 is shifted to zero, only a 2x2 matrix is required
	% rotate and scale,so that x1,y1 -> 0,1
	l2 = P2(:,1).^2 + P2(:,2).^2;
	% sine is 1/l1, but scaling is 1/l1^2
	c  = P2(:,1)./l2;
	s  = P2(:,2)./l2;
	R  = zeros(2,2,n);
	% has to be indeed s and -s
	Ri(1,1,:) =  c;
	Ri(1,2,:) =  s;
	Ri(2,2,:) =  c;
	Ri(2,1,:) = -s;
	P0       = matvec3(Ri,P0.');
%	P2       = matvec3(Ri,P2.');
	P3       = matvec3(Ri,P3.');

	A        = zeros(2,2,n);
	% A(1,1,:) equals P2(1,:) equals 1
	A(1,1,:) = 1;
	A(1,2,:) = P3(1,:);
	% A(2,1,:) equals P2(2,:) equals 0
	A(2,1,:) = 0;
	A(2,2,:) = P3(2,:);
	b        = P0;
	end

	Ai       = inv2x2(A);
	c(2:3,:) = matvec3(Ai,b);

	% shift to origin of X1
	c(1,:) = 1-c(2,:)-c(3,:);
end % tobarycentric2

