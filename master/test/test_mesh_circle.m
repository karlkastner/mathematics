function test_mesh()
	%[p e t]=initmesh('geometry_circle_with_hole','hmax',0.5);
	L0 = 10;
	%[p e t] = initmesh(@(varargin) geometry_rectangle(L0,varargin{:}),'hmax', L0);
	[p e t] = initmesh('geometry_rectangle','hmax', L0);
	pdemesh(p,e,t);
end % test_mesh

function [x y] = grl(varargin)
	L0 = 10;
	[x y] = geometry_rectangle(10,varargin{:})
end

