% Di 26. Jan 20:14:09 CET 2016
%% cartesian to barycentric coordinates
function c = tobarycentric2(P1,P2,P3,P0)
	order = 1;

	% translate so that P1 is at origin for better conditioning
	% TODO scale and roatate to unit triangle

	P0 = P0-P1;
	P2 = P2-P1;
	P3 = P3-P1;
	P1 = zeros(size(P1));

	% set up interpolation matrices
	A    = vander_2d([P1(:,1),P2(:,1),P3(:,1)].', ...
			 [P1(:,2),P2(:,2),P3(:,2)].', order);
	
	b    = vander_2d(P0(:,1).',P0(:,2).',order);

	% TODO exploit that P1=0 and reduce to 2x2 system
	Ai   = inv3x3(A);
	% c := b'*Ai = (Ai'*b)'
	Ait  = transpose3(Ai);
	b    = shiftdim(b,+1);
	c    = matvec3(Ait,b);
	c    = squeeze(c);
end % tobarycentric2

