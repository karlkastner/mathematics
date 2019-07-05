% Fri 12 Jan 16:38:02 CET 2018
%% volume of a tetrahedron
function V = tetra_volume(X,Y,Z)
	A = vander_3d(X',Y',Z',1);
	%A = [X.^0 X Y Z]
	V = 1./6*squeeze(det4x4(A));
end

