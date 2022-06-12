% Wed 11 May 13:10:56 CEST 2022
function [fx, fy, fr, ft, Tx, Ty, mask, N] = fourier_axis_2d(L,n)
	[fx,Tx]=fourier_axis(L(1),n(1));
	[fy,Ty]=fourier_axis(L(2),n(2));
	fr = hypot(cvec(fx),rvec(fy));
	ft = atan2(rvec(fy),cvec(fx));
end

