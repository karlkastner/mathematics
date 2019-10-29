% Wed Aug  8 13:40:23 MSK 2012
% Karl Kästner, Berlin

function plot_mesh(pflag)
	alpha = 1.0;
	k = 6;

	subplot(1,4,1)
	[P T Bc] = mesh_2d_uniform([9 9],[1 1]);
	pdx=1:size(T,1);
	colormap('default');
	c = 0.5 + 0.5*colormap();
%	c(1,:) = c(2,:);
	colormap(c);
%	color = ones(size(pdx))'*c(1,:);
	color = ([2;ones(size(T,1),1)]*[1 1 1]);
	T_ = [T(1,:); T];
	h =fill([ P(T_(:,1),1)  P(T_(:,2),1)  P(T_(:,3),1)]', ...
			[ P(T_(:,1),2)  P(T_(:,2),2)  P(T_(:,3),2)]',color'); %,'facecolor','none','linewidth',2);
	%set(h,'facealpha', alpha);
	set(h,'facecolor','none')
	axis equal
	axis tight
	axis([0 1 0 1.01])
        set(gca,'visible','off');


	
	[P T Bc] = mesh_2d_uniform([3 3],[1 1]);
	%mesh = Mesh_2d(P,T,Bc);
	c = colormap();

	for (idx = 1:k)

	mesh = Mesh_2d(P,T,Bc);
	mesh.element_neighbours();
	N = mesh.N;
	s=sum([ P(T(:,1),:) P(T(:,2),:) P(T(:,3),:) ],2);
	[s fdx] = sort(s);
	M = fdx(1:2);
	[P T Bc] = refine_2d_21(P,T,Bc,N,M);
	pdx=1:size(T,1);
	lw = 1;
	%color = ones(size(T,1),1);
	color = (1:(size(T,1)))';
	subplot(1,4,2)
	h =fill([ P(T(pdx,1),1)  P(T(pdx,2),1)  P(T(pdx,3),1)]', ...
			[ P(T(pdx,1),2)  P(T(pdx,2),2)  P(T(pdx,3),2)]',color'); %,'facecolor','none','linewidth',2);
	%set(h,'facealpha', alpha);
	set(h,'facecolor','none')
	axis equal
	axis tight
	axis([0 1 0 1.01])
        set(gca,'visible','off')
	end

	subplot(1,4,3)
	[P T Bc] = mesh_3d_uniform(9*[1 1 1],1*[1 1 1]);
%	c = 0.5*size(colormap(),1)*ones(size(T,1),1)*[1 2 3];
%	colormap([0 0 1]);
	c = [ones(size(T,1)-1,1);2];
	%h=tetramesh(T,P,c)
	C = (1:size(Bc,1));
	cla();
%		h = patch( [P(Bc(:,1), 1) P(Bc(:,2), 1) P(Bc(:,3), 1) ]', ...
%		       [P(Bc(:,1), 2) P(Bc(:,2), 2) P(Bc(:,3), 2) ]', ...
%		       [P(Bc(:,1), 3) P(Bc(:,2), 3) P(Bc(:,3), 3) ]', C );
	for idx=1:size(Bc,1)
		h = patch( [P(Bc(idx,1), 1) P(Bc(idx,2), 1) P(Bc(idx,3), 1) ]', ...
		       [P(Bc(idx,1), 2) P(Bc(idx,2), 2) P(Bc(idx,3), 2) ]', ...
		       [P(Bc(idx,1), 3) P(Bc(idx,2), 3) P(Bc(idx,3), 3) ]', 1 );
		set(h,'facealpha', alpha);
		set(h,'facecolor','w')
		hold on
	end
        set(gca,'visible','off')

	subplot(1,4,4)
	[P T Bc] = mesh_3d_uniform(3*[1 1 1],1*[1 1 1],0.5*[1 1 1]);
	tree = Tree_3d(P,T,Bc);
	for (idx = 1:k)
		s=sum([ P(T(:,1),:) P(T(:,2),:) P(T(:,3),:) P(T(:,4))],2);
		[s fdx] = sort(s);
		M = fdx(1);
		tree.refine(M);
		mesh = tree.generate_mesh();
		[P T Bc] = get_mesh_arrays(mesh);
	end
	%h=tetramesh(T,P)
%	C = (1:size(Bc,1));
%	h = patch( [P(Bc(:,1), 1) P(Bc(:,2), 1) P(Bc(:,3), 1) ]', ...
%	       [P(Bc(:,1), 2) P(Bc(:,2), 2) P(Bc(:,3), 2) ]', ...
%	       [P(Bc(:,1), 3) P(Bc(:,2), 3) P(Bc(:,3), 3) ]', C );
%	set(h,'facealpha', alpha);
%	set(h,'facecolor','w')

	cla();
	for idx=1:size(Bc,1)
		h = patch( [P(Bc(idx,1), 1) P(Bc(idx,2), 1) P(Bc(idx,3), 1) ]', ...
		       [P(Bc(idx,1), 2) P(Bc(idx,2), 2) P(Bc(idx,3), 2) ]', ...
		       [P(Bc(idx,1), 3) P(Bc(idx,2), 3) P(Bc(idx,3), 3) ]', 1 );
		set(h,'facealpha', alpha);
		set(h,'facecolor','w')
		hold on
	end
        set(gca,'visible','off')

	if (nargin() > 0 && pflag)
		preparePrint();	
		print -depsc ../img/meshes-2d-3d.eps
		system('epstopdf ../img/meshes-2d-3d.eps')
	end

	%colormap (c(10,:))
end
