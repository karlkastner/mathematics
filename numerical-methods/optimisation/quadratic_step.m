% Sa 6. Feb 14:21:33 CET 2016
% Karl Kastner, Berlin
%
%% single step of the quadratic programming
% TODO, line search
function x_opt = quadratic_step(H,g,x0)

[V E] = eig(H);
% check if this a saddle point matrix, if so zero this dimension
me = min(diag(E));
if (me < 0)
	E = diag(E);
	Ei = 1./E;
	Ei(E < 0) = 0;
	Hi = V'*diag(Ei)*V;
	x_opt = x0 - (Hi * g);
	printf('Warning, hessian of objective function is not SPD\n');
	E
else
	% find optimal parameters
	x_opt = x0 - (H \ g);
end

end

