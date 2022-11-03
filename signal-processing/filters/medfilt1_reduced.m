% Mon  9 Jul 09:35:42 CEST 2018
%% median filter with padding
function y = medfilt1_padded(x,n)
	if (isvector(x))
		x = cvec(x);
	end
	y = medfilt1(x,n);
	m = floor(n/2);
	for idx=1:m
		y(idx,:)       = median(x(1:idx*2-1,:));
		y(end-idx+1,:) = median(x(end-2*idx+2:end,:));
	end
end

