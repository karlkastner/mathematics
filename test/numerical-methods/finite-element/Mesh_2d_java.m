% Di 3. Nov 16:26:47 CET 2015
% Karl Kastner, Berlin
% matlab wrapper for the Mesh_2d class

% TODO make refine_2d_21 and get_mesh_arrays part of Mesh_2d
classdef Mesh_2d_java < handle
	properties
		mesh_2d
	end % properties
	methods
	% constructor
	function obj = Mesh_2d_java(P,T,Bc)
		if (nargin() < 2)
			if (strcmp(class(P),'MMesh'))
				mesh = P;
				P  = mesh.point;
				T  = mesh.elem;
				Bc = [mesh.edge(mesh.bnd,1),mesh.edge(mesh.bnd,2)];
			else
				obj.mesh_2d = P;
			end
		else
		% create java object
		obj.mesh_2d = Mesh_2d(P,T,Bc);
		% create neighbourhood relations
		obj.mesh_2d.element_neighbours();
		%obj.mesh.N = Nm;
		end
	end % constructor
	
	function [mmesh obj] = MMesh(obj)
		mmesh         = MMesh();
		mmesh.point   = obj.P; 
		mmesh.elem    = obj.T;  
	end

	function [nt obj] = nt(obj)
		nt = obj.mesh_2d.nt;
	end

	function obj = refine(obj,marked)
		% refinement was not yet ported to java, so fetch matrices and refine in matlab
		[P T Bc Nm] = get_mesh_arrays(obj.mesh_2d);
		[P T Bc Nm] = refine_2d_21(P, T, Bc, Nm, marked);
		% re-create java object
		obj.mesh_2d = Mesh_2d(P,T,Bc);
		% set neighbours
		% TODO, should be accepted by the java constructor
		obj.mesh_2d.N = Nm;
		%obj.mesh_2d.element_neighbours();
	end

	% TODO, search tree (requires parenthood)

	% pseudo members
	% points
	function [P obj] = P(obj)
		P  = double(obj.mesh_2d.P);
		P  = P(1:obj.mesh_2d.np,:);
	end
	% elements (triangles)
	function [T obj] = T(obj)
		T  = double(obj.mesh_2d.T);
		T  = T(1:obj.mesh_2d.nt,:);
	end
	% boundary segments
	function [Bc obj] = Bc(obj)
		Bc = double(obj.mesh_2d.Bc);
		Bc = Bc(1:obj.mesh_2d.nb,:);
	end
	% neighbourhood indicator
	function [N obj] = N(obj)
		N  = double(obj.mesh_2d.N);
	end
	end % methods
end % class Mesh_2d_java

