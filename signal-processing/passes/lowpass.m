% 2016-02-27 16:37:59.183617250 +0100
%% low pass filter
% TODO, bad, bad, bad this is not an ideal cos-filter
% f = f_c/f_s = dt/T_c
function f = lowpass(n,order) %f_c,order)
%	n = sqrt(4*0.196/f_c^2 + 1);
%	n = round(n);
	f0 = ones(n,1);
	f  = f0;
	for idx=1:order-1
		f = conv(f,f0);
	end
	f = f/sum(f);
end

