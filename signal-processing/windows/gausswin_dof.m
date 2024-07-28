% 2024-07-27 15:20:59.131050614 +0200
% effective sample size of a gaussian window on an uncorrelated dataset
function dof gaussian_window_dof(s,dx)
	dof = 2*sqrt(pi)*s/dx;
end

