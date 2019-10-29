% Thu Aug  9 21:15:11 MSK 2012
% Karl Kästner, Berlin

function plot_erro_estimation(pflag)

	lw = 1;
	ms = 16;
	z = 1.7;
	alpha=0.25;
	colour = [0.5 0.5 0.5];

	figure(1);
	clf();
	
	P_ = [ 0 0; 1 0; 0.5 sqrt(1-0.25) ];
	T_ = [ 1 2 3];
	Bc_ = [2 3; 1 3; 1 2];

	mesh = Mesh_2d(P_,T_,Bc_);
	mesh.promote(2);
	mesh.demote();
	mesh_ = Mesh_2d(mesh.P,mesh.S,mesh.Bc);
	
	subplot(1,4,1)
	P = mesh.P; T=mesh.T; Bc = mesh.Bc;
	h =fill([ P(T(:,1),1)  P(T(:,2),1)  P(T(:,3),1)]', ...
			[ P(T(:,1),2)  P(T(:,2),2)  P(T(:,3),2)]',1,'linewidth',lw);
	set(h,'facecolor','none');
	hold on;
	plot(P(:,1), P(:,2),'k.','Markersize',ms);
	set(gca,'visible','off')
	axis equal
	
	subplot(1,4,2)
	P = mesh_.P; T=mesh_.T; Bc = mesh_.Bc;
	h =fill([ P(T(:,1),1)  P(T(:,2),1)  P(T(:,3),1)]', ...
			[ P(T(:,1),2)  P(T(:,2),2)  P(T(:,3),2)]',1,'linewidth',lw);
	set(h,'facecolor','none');
	hold on;
	plot(P(:,1), P(:,2),'k.','Markersize',ms);
	set(gca,'visible','off')
	axis equal


	mesh = Mesh_2d(P_,T_,Bc_);
	mesh.promote(4);
	mesh.demote();
	mesh_ = Mesh_2d(mesh.P,mesh.S,mesh.Bc);
	
	subplot(1,4,3)
	P = mesh.P; T=mesh.T; Bc = mesh.Bc;
	h =fill([ P(T(:,1),1)  P(T(:,2),1)  P(T(:,3),1)]', ...
			[ P(T(:,1),2)  P(T(:,2),2)  P(T(:,3),2)]',1,'linewidth',lw);
	set(h,'facecolor','none');
	hold on;
	plot(P(:,1), P(:,2),'k.','Markersize',ms);
	set(gca,'visible','off')
	axis equal
			
	subplot(1,4,4)
	P = mesh_.P; T=mesh_.T; Bc = mesh_.Bc;
	h =fill([ P(T(:,1),1)  P(T(:,2),1)  P(T(:,3),1)]', ...
			[ P(T(:,1),2)  P(T(:,2),2)  P(T(:,3),2)]',1,'linewidth',lw);
	set(h,'facecolor','none');
	hold on;
	plot(P(:,1), P(:,2),'k.','Markersize',ms);
	set(gca,'visible','off')
	axis equal

	figure(2);
	clf();	

	P = [0 0 0;
	     1 0 0;
	     0 1 0;
	     0 0 1];
	P_ = [0 0 0;
             1 0 0;
             0.5 sqrt(1-0.5^2) 0;
             0.5 0.5*sqrt(1-0.5^2) sqrt(1 - 0.25 -0.25*(1-0.25))];
	P_ = P_ - ones(4,1)*mean(P_);
	T_ = [1 2 3 4];
	Bc_ = [2 3 4;
	       1 3 4;
	       1 2 4;
	       1 2 3];

	mesh = Mesh_3d(P_,T_,Bc_);
	mesh.promote(2);
	mesh.demote();

	subplot(1,4,1)
	[P T Bc] = get_mesh_arrays(mesh);
	h = tetramesh(T, P, 1,'linewidth',lw);
	set(h,'facecolor',colour);
	set(h,'facealpha',alpha)
	set(gca,'visible','off')
	hold on;
	plot3(P(:,1), P(:,2), P(:,3), 'k.','Markersize',ms);
	zoom(z)

	subplot(1,4,2)
	[P T Bc] = get_mesh_arrays(mesh); T=mesh.S;
	h = tetramesh(T, P, ones(8,1),'linewidth',lw);
	set(h,'facecolor',colour);
	set(h,'facealpha',alpha)
	set(gca,'visible','off')
	hold on;
	plot3(P(:,1), P(:,2), P(:,3), 'k.','Markersize',ms);
	zoom(z)

	mesh = Mesh_3d(P_,T_,Bc_);
	mesh.promote(4);
	mesh.demote();

	subplot(1,4,3)
	[P T Bc] = get_mesh_arrays(mesh);
	h = tetramesh(T, P, 1,'linewidth',lw);
	set(h,'facecolor',colour);
	set(h,'facealpha',alpha)
	set(gca,'visible','off')
	hold on;
	plot3(P(:,1), P(:,2), P(:,3), 'k.','Markersize',ms);
	zoom(z)

	subplot(1,4,4)
	[P T_ Bc] = get_mesh_arrays(mesh); T=mesh.S;
	h = tetramesh(T, P, ones(size(T,1),1),'linewidth',lw);
	set(h,'facecolor',colour);
	set(h,'facealpha',alpha)
	set(gca,'visible','off')
	hold on;
	plot3(P(:,1), P(:,2), P(:,3), 'k.','Markersize',ms);
	zoom(z)

	if (nargin() > 0 && pflag)
		figure(1);
		print -deps ../img/error-estimation-2d.eps
		system('epstopdf ../img/error-estimation-2d.eps')
		figure(2);
		print -deps ../img/error-estimation-3d.eps
		system('epstopdf ../img/error-estimation-3d.eps')
	end


%{	
	figure(3); clf
	for idx=1:size(T,1)
		clf
		h = tetramesh(T(idx,:), P, 1,'linewidth',lw); hold on
		pdx=T(idx,:);
		plot3(P(pdx,1), P(pdx,2), P(pdx,3), 'k.','Markersize',ms);
		set(h,'facecolor','none');
		pause(1)
	end
%}
%	for idx=1:size(P,1)
%		text(P(idx,1), P(idx,2), P(idx,3), num2str(sdx(idx)-1),'FontSize',24); hold on
%	end
%
%		for jdx=1:size(P,1)
%			line([0 P(jdx,1)], [0 0], [0 0])
%			line([P(jdx,1) P(jdx,1)], [0 P(jdx,2)], [0 0])
%			line([P(jdx,1) P(jdx,1)], [P(jdx,2) P(jdx,2)], [0 P(jdx,3)])
%		end

%[T sdx] = sort(T)
%A=reshape(sdx(Bc(:)),4,15)-1
%size(unique(A(:)))



end

