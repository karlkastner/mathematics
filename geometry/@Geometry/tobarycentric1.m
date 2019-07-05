% Sun  5 Feb 14:49:56 CET 2017
%% cartesian to barycentric coordinates
function c = tobarycentric1(P1,P2,P0)
	order = 1;

	% shift to local coordinates for better conditioning
	P0 = P0-P1;
	P2 = P2-P1;
	P1 = 0.*P1;

	% scale to unity
	% P0=0/L;
	L  = P2; % == P2-P1==P2-0

	P0 = P0./L;
	P2 = P2./L;

	% TODO eX is now [0,1] exploit this
	%A  = vander_1d(eX',order);
	% set up interpolation matrices
	A  = vander_1d([P1,P2]', order);
	% this is now [1 0; 1 1]
	A0   = vander_1d(P0', order)';
	Ai   = inv2x2(A);
	% C := A'*Ai = (Ai'*A)'
	Ait  = transpose3(Ai);
	%Ait  = transpose2x2(Ai);
	c    = matvec3(Ait,A0);
	c    = squeeze(c);
end % tobarycentric1

