% Fri Jul 27 17:07:53 MSK 2012
% Karl Kästner, Berlin
%
% TODO exact integration with f=1
%
% T : global point indices (size = [1 6])
% P : point coordinates (size)
% f : potential function (function handle)
function [A buf] = assemble_3d_phi_phi_java(mesh, func, int)

	% load the integration sample point coordinates for standard triangle
	[w b flag] = feval(int);

	% get matrix entries
		%asm = Assemble_2d_phi_phi;
		%buf = asm.assemble(mesh, func, w, b, flag);
	buf = mesh.assemble_phi_phi(func, w, b, flag);

	% construct the matrix from the buffer
	A = sparse(buf(:,1), buf(:,2), buf(:,3), double(mesh.np), double(mesh.np));
	% complete symmetry
	A = A + A';

end % assemble_3d_phi_phi()

