% Sun Jul 29 17:09:59 MSK 2012
% Karl Kästner, Berlin

function plot_basis_functions(pflag)
	clf()
	lw=1;
	ms = 12;
	colormap([0.06125 0.25 1])
	cols = 6;
	rows = 4;

	% lines

	for idx=1:5
		subplot(2*rows,cols,cols+idx);
		plot(2*(0:(idx))/(idx+1), zeros(6,(idx+1)), 'k.-', 'linewidth', lw, 'Markersize', ms);
		text(-0.0,-0.5,['p =' num2str(idx) ' n = ' num2str((idx+1))])
	        set(gca,'visible','off')
	%	axis equal;
		axis tight
	end

	% triangles

	T_ = [1 2 3];
	P_ = [0 0
	     0 1;
	     1 0];
	Bc_ = [1 2 -1;
	      2 3 -1;
	      1 3 -1];
%	[P T Bc] = mesh_2d_uniform([2 2], [1 1], [0 0]);
	
	alpha = 0.4;
	for idx=1:5
		mesh=Mesh_2d(P_,T_,Bc_);
		mesh.promote(idx);
		P = mesh.P;
		T = T_;

		subplot(rows,cols,cols+idx)
		cla
		pdx=1;
	%	h =patch([ P(T(pdx,1),1)  P(T(pdx,2),1)  P(T(pdx,3),1)]', ...
	%			[ P(T(pdx,1),2)  P(T(pdx,2),2)  P(T(pdx,3),2)]',[0.5 1 0.5]);
		h =fill([ P(T(pdx,1),1)  P(T(pdx,2),1)  P(T(pdx,3),1)]', ...
				[ P(T(pdx,1),2)  P(T(pdx,2),2)  P(T(pdx,3),2)]',1,'linewidth',lw);
	%	,'FaceColor','g'); hold on
		hold on
		set(h,'facealpha', alpha);
		set(h,'facecolor','none');
	%	h = patch(T(1,:), P, 'r');
		plot(P(:,1),P(:,2),'.k','Markersize',ms)
		    set(gca,'visible','off')
		text(-0.0,-0.1,['p =' num2str(idx) ' n = ' num2str((idx+1)*(idx+2)/2)])
	%	    text(-0.6,-0.2,num2str((idx+1)*(idx+2)/2));
		a = axis;
		da = [-0.1 0.1 -0.1 0.1];
		axis(a+da)
		axis equal;
		axis tight
	
	end

	alpha = 0.25;
	T_ = [1 2 3 4];
	P_ = [0 0 0;
	     1 0 0;
	     0 1 0;
	     0 0 1];
	Bc_ = [1 2 3;
	      1 2 4;
	      1 3 4;
	      2 3 4];
	L0 = 1;
	n = 2;
	%[P_ T_ Bc_] = mesh_3d_uniform([n n n], [L0 L0 L0]);
%	P_ = [0 0 0;
%             1 0 0;
%             0.5 sqrt(1-0.5^2) 0;
%             0.5 0.5*sqrt(1-0.5^2) sqrt(1 - 0.25 -0.25*(1-0.25))];

	colormap([0.06125 0.25 1])
	colour = [0.5 0.5 0.5];
	colour = 'none';
	alpha = 0.5;
	z = 1;
	sp = 0; %-0.05;

	for idx=1:5
		mesh = Mesh_3d(P_,T_,Bc_);
		mesh.promote(idx);
		[P T Bc] = get_mesh_arrays(mesh);
		subplot(rows,cols,2*cols+idx)
		h = tetramesh(T(1,:), P, 1,'linewidth',lw);
		hold on
		set(h,'facecolor', colour);
		set(h,'facealpha', alpha);
	        set(gca,'visible','off');
		plot3(P(:,1),P(:,2),P(:,3),'.k','MarkerSize',ms);
		text(1.75, 0.25, ['p =' num2str(idx) ' n = ' num2str((idx+1)*(idx+2)*(idx+3)/6)]);
		%text(-0.5, 1, ['p =' num2str(idx) ' n = ' num2str((idx+1)*(idx+2)*(idx+3)/6)]);
		%view(36,-22.5);
		view(102.5,12.5)
		axis equal;
		axis tight;
		zoom(z);
		hold off;

%		for jdx=1:size(P,1)
%			line([0 P(jdx,1)], [0 0], [0 0])
%			line([P(jdx,1) P(jdx,1)], [0 P(jdx,2)], [0 0])
%			line([P(jdx,1) P(jdx,1)], [P(jdx,2) P(jdx,2)], [0 P(jdx,3)])
%		end


	end
	for idx=1:5
		subplot(rows,cols,3*cols-idx);
		pos = get(gca,'Position');
		pos(1) = pos(1)+sp;
		set(gca,'Position',pos);
	end


	if (nargin() > 0 && pflag)
		preparePrint();
		print -depsc ../img/dof-2d-3d.eps
		system('epstopdf ../img/dof-2d-3d.eps')
	end
end

