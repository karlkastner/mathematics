% Tue May  1 01:19:46 MSK 2012
% Karl Kästner, Berlin
%
% integrates product of the derivatives of the test functions
%
% T : points on the triangle point, first three points are corners (size = [1 .. 6])
function [A buf] = assemble_2d_dphi_dphi_java(mesh, func, int)

	% load the integration sample point coordinates for standard triangle
	[w b flag] = feval(int);

	% get matrix entries	
	asm = javaObject('Assemble_2d_dphi_dphi');
	buf = asm.assemble(mesh, func, w, b);
	if (exist('OCTAVE_VERSION') > 0)
		% TODO FIXME
		buf = java2mat(buf);
	end
	% construct matrix from buffer
	A = sparse(buf(:,1), buf(:,2), buf(:,3), size(mesh.P,1), size(mesh.P,1));
	% complete symmetry
	A = A + A';
	
end % function assemble_2d_dphi_dphi

