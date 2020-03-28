% Sun Aug  5 23:38:09 MSK 2012
% Karl Kästner, Berlin

function [P T Bc N] = get_mesh_arrays(mesh)
	P  = double(mesh.P);
	P  = P(1:mesh.np,:);
	T  = double(mesh.T);
	T  = T(1:mesh.nt,:);
	Bc = double(mesh.Bc);
	Bc = Bc(1:mesh.nb,:);
	N  = double(mesh.N);
% do not sort, otherwise indices on tree not same
%	Bc = sortrows(sort(Bc,2));
%	T  = sortrows(sort(T,2));
end

