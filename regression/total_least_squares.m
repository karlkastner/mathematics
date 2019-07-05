% Thu Aug 22 07:20:43 UTC 2013
% Karl Kästner, Berlin
%
%% total least squares
%
function [B, E, F] = tls(X,Y)
	% (X + E)B = (Y+F), min ||[E F]||_F
	% [X Y] = [Ux Uy] * [Sx  0; 0  Sy] * [Vxx Vxy; Vyx Vyy]'; 
	nx = size(X,2);
	[U S V] = svd([X Y],0);
	S = diag(S);
	Uy  = U(:,nx+1:end);
	Sy  = diag(S(nx+1:end));
	Vxy = V(1:nx,nx+1:end);
	Vyy = V(nx+1:end,nx+1:end);
	EF = -Uy * Sy * [Vxy; Vyy]';
	E = EF(:,1:nx);
	F = EF(:,nx+1:end);
	B = -Vxy / Vyy; % rdevide is correct
end

