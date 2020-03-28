	n = 20*[1 1];
	L0 = 1*[1 1];
	%x0 = [0.25 1/3];
	x0 = 0.5*L0;

% todo, x0 into mesh_2d_uniform

	[P T Bc X] = mesh_2d_uniform(n+1, L0);
	P(:,1) = P(:,1) - x0(1);
	P(:,2) = P(:,2) - x0(2);

	%K_r = 0;
	%K_t = 0.1;
	K_t = [0 0];
	clf
	flag = 1;
	subplot(2,2,1)
	K_r = [1 0 0];
	K_t = [0 0];
	P_ = distort(P,K_r,K_t,L0,flag,2);
	[area l_boundary h_side s_angle C] = regularity_2d(P_,T,Bc);
	[min(min(h_side)) min(min(s_angle)) min(min(h_side./s_angle))]
	%P_ = P_/sqrt(2);
	display_2d(P_,T,Bc,0,[],[]);
	axis equal, axis tight

	subplot(2,2,2)
	K_r = [1 1 0];
	K_t = [0 0];
	P_ = distort(P,K_r,K_t,L0,flag,4);
	[area l_boundary h_side s_angle C] = regularity_2d(P_,T,Bc);
	[min(min(h_side)) min(min(s_angle)) min(min(h_side./s_angle))]
	%[min(min(h_side)) mean(mean(h_side)) min(min(s_angle))]
	display_2d(P_,T,Bc,0,[],[]);
	axis equal, axis tight

	[P T Bc X] = mesh_2d_uniform(2*n+1, L0);
	P(:,1) = P(:,1) - x0(1);
	P(:,2) = P(:,2) - x0(2);

	subplot(2,2,3)
	K_r = [0 0 1];
	K_t = [0 0];
	P_ = distort(P,K_r,K_t,L0,flag,2);
	[area l_boundary h_side s_angle C] = regularity_2d(P_,T,Bc);
	[min(min(h_side)) min(min(s_angle)) min(min(h_side./s_angle))]
%	[min(min(h_side)) mean(mean(h_side)) min(min(s_angle))]
	display_2d(P_,T,Bc,0,[],[]);
	axis equal, axis tight



	subplot(2,2,4)
	K_r = [1 0 0];
	K_t = [0 0];
	P_ = distort(P,K_r,K_t,L0,flag,4);
	[area l_boundary h_side s_angle C] = regularity_2d(P_,T,Bc);
	[min(min(h_side)) min(min(s_angle)) min(min(h_side./s_angle))]
%	[min(min(h_side)) mean(mean(h_side)) min(min(s_angle))]
	display_2d(P_,T,Bc,0,[],[]);
	axis equal, axis tight
%{
	P
	P = P + 0.1*rand(size(P));
	
	for bdx=1:size(B,1)
		p1 = B(bdx,1);
		boundary =
		p1 = B(bdx,1);
		P(
function P = project_rectangle(pdx,P,bdx,L0,x0)
	P = project_rectangle	
%}

