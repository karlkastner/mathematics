% Sat 14 Jan 15:02:43 CET 2023
function [Rt,ti] = autocorr_angular(R,L,m)
	n = size(R);
	[x,y,r,t] = fourier_axis_2d(L,n);
	Lmax = min(L)/2;
	ti = innerspace(-pi,pi,m);
	dt = ti(2)-ti(1);
	id = floor((t + pi)/dt)+1;
	fdx = abs(r)<Lmax;
	Rt = accumarray(id(fdx),R(fdx),[length(ti),1],@mean, 0);
end

