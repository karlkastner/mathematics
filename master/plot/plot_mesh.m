function plot_mesh()
	javaaddpath('.')
	[P T Bc] = mesh_3d_uniform(-1,1);
	P = [0 0 0;
             1 0 0;
             0.5 sqrt(1-0.5^2) 0;
             0.5 0.5*sqrt(1-0.5^2) sqrt(1 - 0.25 -0.25*(1-0.25))];
	P = P - ones(4,1)*mean(P);
	P
	tetramesh(T,P)

	tree = Tree_3d(P,T,Bc);
	tree.refine(1);
	mesh = tree.generate_mesh();
	P = mesh.P; P = P(1:mesh.np,:);
	T = mesh.T; T = T(1:mesh.nt,:);
	B = mesh.Bc; B = B(1:mesh.nb,:);


	subplot(1,2,2)
	[P T] = explode(P,T,[],1.25);
P
T
	colormap('default')
	c = colormap
	p = 0.06125;
	colormap((1-p)*c + p*ones(size(c,1),1)*[ 0 0 1]);
	h=tetramesh(T,P)
	alpha = 0.5;
	set(h,'facealpha', alpha);
%	view(60,30);
%	view(90,30);
	view(80,10)
        set(gca,'visible','off')
	axis equal, axis tight

	subplot(1,2,1)
	[P T Bc] = mesh_2d_uniform([2 2],[1 1]);
	%mesh = Mesh_2d(P,T,Bc);
	N = FEM.elementNeighbours(P,T,Bc);
	[P T Bc] = refine_2d_21(P,T,Bc,N,1)
	pdx=1:size(T,1);
	lw = 1;
	color = [1 2 1 1 1 2]'
	h =fill([ P(T(pdx,1),1)  P(T(pdx,2),1)  P(T(pdx,3),1)]', ...
			[ P(T(pdx,1),2)  P(T(pdx,2),2)  P(T(pdx,3),2)]',color','facecolor','none','linewidth',2);
	set(h,'facealpha', alpha);
	axis equal
	axis tight
        set(gca,'visible','off')

	preparePrint();	
	print -depsc ../../img/tetra-tesselation.eps
	system('epstopdf ../../img/tetra-tesselation.eps')
	
end

% todo, mesh plots
% -todo : refinement plots

