% Tue May  1 01:20:41 MSK 2012
% Karl Kästner, Berlin
%
% TODO exact integration with f=1
%
% T : global point indices (size = [1 6])
% P : point coordinates (size)
% f : potential function (function handle)
function [A buf] = assemble_2d_phi_phi_java(mesh, func, int)

	% load the integration sample point coordinates for standard triangle
	[w b flag] = feval(int);

	% get matrix entries
	asm = javaObject('Assemble_2d_phi_phi');
	buf = asm.assemble(mesh, func, w, b, flag);

	if (exist ('OCTAVE_VERSION') > 0)
		% TODO FIXME
		buf = java2mat(buf);
	end

	% construct the matrix from the buffer
	A = sparse(buf(:,1), buf(:,2), buf(:,3), size(mesh.P,1), size(mesh.P,1));
	% complete symmetry
	A = A + A';

end % assemble_2d_phi_phi

