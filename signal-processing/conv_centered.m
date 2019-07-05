% Wed 18 Jan 09:51:28 CET 2017
% Karl Kastner, Berlin
%% convolve x with filter window f
%% when length of f is even, this guarantees a symmetric result (no off by on
%% displacement) by making the lenght of f odd at first
function xf = conv_centered(x,f)
	flag = 0;
	f  = window_make_odd(f,5);
	ir = isrow(x);
	if (ir)
		x = x.';
	end
	
	xf = zeros(size(x));	

	for idx=1:size(x,2)
		xf(:,idx)= conv(x(:,idx),f,'same');
	end

	if (ir)
		xf = xf.';
	end
end

