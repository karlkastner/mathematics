% Wed Aug  8 14:43:49 MSK 2012
% Karl Kästner, Berlin

function plot_error_estimation()
	clf();
	lw = 1;
	ms = 12;
        % p=1, refined
	subplot(2,4,1)
	[P T Bc] = mesh_2d_uniform(-1,[]);
	mesh = Mesh_2d(P,T,Bc);
	mesh.element_neighbours();
	N = mesh.N;
	[P T Bc] = refine_2d_21(P,T,Bc,N,1);
	h =fill([ P(T(:,1),1)  P(T(:,2),1)  P(T(:,3),1)]', ...
			[ P(T(:,1),2)  P(T(:,2),2)  P(T(:,3),2)]',1,'linewidth',lw);
	hold on
	set(h,'facecolor','none');
		plot(P(:,1),P(:,2),'.k','Markersize',ms)

	% p=2
	subplot(2,4,2)
	[P T Bc] = mesh_2d_uniform(-1);
	mesh = Mesh_2d(P,T,Bc);
	mesh.element_neighbours();
	mesh.prefetch();
	mesh.promote_3_6();
	P = mesh.P;
	T = mesh.T;
	Bc = mesh.Bc;
	N = mesh.N;
	h =fill([ P(T(:,1),1)  P(T(:,2),1)  P(T(:,3),1)]', ...
			[ P(T(:,1),2)  P(T(:,2),2)  P(T(:,3),2)]',1,'linewidth',lw);
	set(h,'facecolor','none');
	hold on
		plot(P(:,1),P(:,2),'.k','Markersize',ms)

        % p=1, refined
	subplot(2,4,1)
	[P T Bc] = mesh_2d_uniform(-1,[]);
	mesh = Mesh_2d(P,T,Bc);
	mesh.element_neighbours();
	mesh.promote_3_10();
	N = mesh.N;
	[P T Bc] = refine_2d_21(P,T,Bc,N,1);
	h =fill([ P(T(:,1),1)  P(T(:,2),1)  P(T(:,3),1)]', ...
			[ P(T(:,1),2)  P(T(:,2),2)  P(T(:,3),2)]',1,'linewidth',lw);
	hold on
	set(h,'facecolor','none');
		plot(P(:,1),P(:,2),'.k','Markersize',ms)

	% p=2
	subplot(2,4,2)
	[P T Bc] = mesh_2d_uniform(-1);
	mesh = Mesh_2d(P,T,Bc);
	mesh.element_neighbours();
	mesh.prefetch();
	mesh.promote_3_6();
3 6 10 15 21 28
	mesh.promote_3_10();
	P = mesh.P;
	T = mesh.T;
	Bc = mesh.Bc;
	N = mesh.N;
	h =fill([ P(T(:,1),1)  P(T(:,2),1)  P(T(:,3),1)]', ...
			[ P(T(:,1),2)  P(T(:,2),2)  P(T(:,3),2)]',1,'linewidth',lw);
	set(h,'facecolor','none');
	hold on
		plot(P(:,1),P(:,2),'.k','Markersize',ms)
	
end

