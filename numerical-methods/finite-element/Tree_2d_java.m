% Di 9. Feb 16:48:58 CET 2016
% 
classdef Tree_2d_java < handle
	properties
		tree_2d
	end
	methods
	function obj = Tree_2d_java(mmesh)
%		mmesh.point = [X Y];
%		mmesh.elem  = elem;
		mmesh.edges_from_elements();	
		Bc            = mmesh.edge(mmesh.bnd,:);
		obj.tree_2d   = Tree_2d(mmesh.point(:,1:2),mmesh.elem,Bc);
	end
	function obj = refine(obj,M)
		obj.tree_2d.fill_R();
		obj.tree_2d.refine(M);
	end
	function [mmesh obj] = generate_MMesh(obj)
		mesh_2d_java  = obj.generate_Mesh_2d();
		mmesh          = mesh_2d_java.MMesh;
	end
	function [mesh_2d_java obj] =generate_Mesh_2d(obj)
		mesh_2d       = obj.tree_2d.generate_mesh(true);
		mesh_2d_java  = Mesh_2d_java(mesh_2d);
	end
	end
end

