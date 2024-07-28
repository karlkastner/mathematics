% 2024-07-27 18:04:21.852429159 +0200
function sd = gausswin_dof2std(dof,dx)
	sd = dof*dx/(2*sqrt(pi));
end

