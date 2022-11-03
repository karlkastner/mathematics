function [x] = bandpass1d(x,p,k)
	if (length(p) == 1)
		p =[p,p];
	end
	if (nargin()<3)
		k = 1;
	end
	for idx=1:k
		%x = meanfilt1(x,5);
		%x = x-meanfilt1(x,20);
		x = filter(1-p(1),[1,-p(1)],x);
		x = flipud(filter(1-p(1),[1,-p(1)],flipud(x)));
		x=  x-filter(1-p(2),[1,-p(2)],x);
		x=  x-flipud(filter(1-p(2),[1,-p(2)],flipud(x)));
	end
end

