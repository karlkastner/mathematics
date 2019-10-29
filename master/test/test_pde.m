function test_pde()

[p,e,t]=initmesh('lshapeg','Hmax',0.2); 
[p,e,t]=refinemesh('lshapeg',p,e,t); 
u=assempde('lshapeb',p,e,t,1,0,1); 
pdesurf(p,t,u)

	a = 1;
	c = @cfunc;
	f = @ffunc;
	a = @afunc;
	a = 1;
	f = 1;
	%c = 1;

	%[ar,g1x,g1y,g2x,g2y,g3x,g3y] = pdetrg(p,t)
	[p e t] = initmesh('geometry_rectangle','hmax', 2);
	[K,M,F,Q,G,H,R] = assempde('boundary_rectangle',p,e,t,c,a,f);

%	[p e t] = initmesh('geometry_circle_with_hole','hmax', 10);
%	u=assempde('boundary_cirlce',p,e,t,c,a,f)
%	[K,M,F,Q,G,H,R] = assempde('boundary_circle',p,e,t,c,a,f);
	full(K)
	full(M)
	eig(full(K),full(M))
	eigs(K,M,1,'SM')
end

function a = afunc(x,varargin)	
	a = ones(1,size(x,2));
end

function f = ffunc(x,varargin)	
	f = ones(1,size(x,2));
end

function c = cfunc(x,varargin)
	c = 1./sqrt(sum(x.^2,1))
%	c = diag(c);
%	c = c';
%	c = [c;c]
end

