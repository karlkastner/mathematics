% 2016-02-27 12:55:45.388567143 +0100
%% bandpass filter
function f = bandpass(m,n,order)
	% filter length is (2*o*(n-1)+1)
	f1 = lowpass(m,order);
	f2 = highpass(n,order);
	f = conv(f1,f2);
end

