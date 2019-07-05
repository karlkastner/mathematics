% Thu 21 Jul 11:18:27 CEST 2016
% Karl Kastner, Berlin
%
%% expand values of fourier series
%
function y = fourier_expand(a,b,x)
	y = a(1)*ones(size(x));
	for idx=1:length(a)-1
		y = y + a(idx+1)*cos(idx*x) + b(idx+1)*sin(idx*x);
	end
end

