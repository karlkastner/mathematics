% 2016-02-27 12:35:36.077558230 +0100 
%% high pass filter
function f = highpass(n,order)
	f = lowpass(n,order);
%	m = (n+1)/2;
	m = (length(f)+1)/2;
	f(m) = 0;
	f = -f/sum(f);
	f(m) = 1;
%	f0 = f;
%	for idx=1:order-1
%		f = conv(f,f0);
%	end
%	f = f/sum(f);
	
end

