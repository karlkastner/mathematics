% 2016-02-27 14:19:53.767431036 +0100
%% frequncy response of a boxcar filter
function h = freqz_boxcar(f,n)
	h = 2/n*abs(sin(2*pi*f*n/2)./sin(2*pi.*f));
end
