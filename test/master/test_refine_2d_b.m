javaaddpath('.');

%{
P = [ 0 0;
      1 0;
      0 1;
      0 -1];
T = [1 2 3;
     2 1 4];
Bc = [1 3 1;
      2 3 1;
      1 4 1;
      2 4 1];
Nm = double(FEM.elementNeighbours(P,T,Bc));
G = zeros(size(T,1),1);

%clf
subplot(2,2,1);
mesh = Mesh(P,T,Bc);
display_2d(mesh, 7, [], [], 'EdgeColor', 'k');
axis equal
%[P T Bc Nm G] = refine_2d(P, T, Bc, 1, Nm, G);
%[P T Bc Nm] = refine_2d_21(P, T, Bc, Nm, 1);
Nm_ = double(FEM.elementNeighbours(P,T,Bc));
norm(Nm-Nm_)
[(1:size(T,1))' Nm Nm_]
subplot(2,2,2);
mesh = Mesh(P,T,Bc);
display_2d(mesh, 7, [], [], 'EdgeColor', 'k');
axis equal
%[P T Bc Nm G] = refine_2d(P, T, Bc, [6;3], Nm, G);
%[P T Bc Nm] = refine_2d_21(P, T, Bc, Nm, 1);
% problem 3,4

subplot(2,2,3);
mesh = Mesh(P,T,Bc);
display_2d(mesh, 7, [], [], 'EdgeColor', 'k');
axis equal
%[P T Bc N G] = refine_2d(P, T, Bc, 1, N, G);
%}
%clf
L0=1;
n=3;
[P T Bc] = mesh_2d_uniform(n*[1 1], 2*L0*[1 1]);
Nm = double(FEM.elementNeighbours(P,T,Bc));
	subplot(2,2,4);
	mesh = Mesh(P,T,Bc);
	display_2d(mesh, 7, [], [], 'EdgeColor', 'k');
T
tic
P=P-L0/2;
Tr = [];
N_ = [];
figure(1)
for idx=1:11
%	pdx = [Bc(:,1); Bc(:,2)];
%	P(pdx,:) = project_circle(P(pdx,:),L0);
%	lm = randi(ceil(size(T,1)/4));
%	M = unique(randi(size(T,1),lm,1)');
	%M = (1:size(T,1));
	%M = [1 4 6 8 3 5];
	%M = 1:8;
	%M = 1:size(T,1);
	M=1;
	subplot(2,2,3); cla;
	mesh = Mesh(P,T,Bc);
	display_2d(mesh, 7, [], [], 'EdgeColor', 'k'); axis equal
	tic();
	[P T Bc Nm] = refine_2d_21(P, T, Bc, Nm, M);
	Tr(idx) = toc()
	N_(idx) = size(P,1);

	mesh = Mesh(P,T,Bc);
	subplot(2,2,4);	cla; display_2d(mesh, 7, [], [], 'EdgeColor', 'k'); axis equal
	Nm_ = double(FEM.elementNeighbours(P,T,Bc));
	if (norm(Nm-Nm_) > 0)
		[(1:size(T,1))' Nm Nm_ Nm-Nm_]
		'error'
		return
	end
	
end
toc
size(P,1)
figure(2);
loglog(N_,Tr);


