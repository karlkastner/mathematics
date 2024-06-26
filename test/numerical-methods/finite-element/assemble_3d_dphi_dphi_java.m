% Fri Jul 27 17:06:18 MSK 2012
% Karl Kästner, Berlin
%
% integrates product of the derivatives of the test functions
%
% T : points on the triangle point, first three points are corners (size = [1 .. 6])
function [A buf] = assemble_3d_dphi_dphi_java(mesh, func, int)
	% load the integration sample point coordinates for standard triangle
	[w b flag] = feval(int);

	% get matrix entries	
		%asm = Assemble_3d_dphi_dphi;
		%buf = asm.assemble(mesh, func, w, b);
	buf = mesh.assemble_dphi_dphi(func, w, b);
	% construct matrix from buffer
	A = sparse(buf(:,1), buf(:,2), buf(:,3), mesh.np, mesh.np);
	% complete symmetry
	A = A + A';
	
end % function assemble_3d_dphi_dphi

