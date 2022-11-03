% Mon  9 Jul 09:35:42 CEST 2018
%% median filter with padding
function x = medfilt1_padded(x,n)
	if (isvector(x))
		x = cvec(x);
	end
	x = [ x(1:n,:); x; x(end-n+1:end,:) ];
	x = medfilt1(x,n);
	x = x(n+1:end-n,:);
end
