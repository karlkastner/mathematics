% Thu 11 Jan 14:55:09 CET 2018
%% cartesian to barycentric coordinates
function c = tobarycentric3(P1,P2,P3,P4,P0)
	order = 1;

	% shift to local coordinates for better conditioning
	% TODO scale and roatate to unit tetra
	P0   = P0-P1;
	P2   = P2-P1;
	P3   = P3-P1;
	P4   = P4-P1;
	P1   = zeros(size(P1));

	% set up interpolation matrices
	A   = vander_3d([P1(:,1),P2(:,1),P3(:,1),P4(:,1)].', ...
			[P1(:,2),P2(:,2),P3(:,2),P4(:,2)].', ...
			[P1(:,3),P2(:,3),P3(:,3),P4(:,3)].', order);
	b   = vander_3d(P0(:,1).',P0(:,2).',P0(:,3).', order);

	% TODO, this can be simplified to 3x3, as P0 has been shifted to 0
	Ai   = inv4x4(A);
	% C := A'*Ai = (Ai'*A)'
	Ait  = transpose3(Ai);
	b    = shiftdim(b,+1);
	c    = matvec3(Ait,b);
	c    = squeeze(c);
end % tobarycentric3

