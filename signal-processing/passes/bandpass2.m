% 2016-02-27 15:46:30.584686222 +0100
% Karl Kastner, Berlin
%% bandpass filter
function f = bandpass2(n,f1,f2,fs)
	f = sort([f1 f2]);
	scale = 1+sqrt(eps);
	if (f2/f1 < scale)
		warning('Freuqencies too close, increasing notch width to machine precission');
		f1 = f1/scale;
		f2 = f2*scale;
	end
	% normalise sample frequencies
	f1 = 2*f1/fs;
	f2 = 2*f2/fs;
	% design filter
	f = fir1(n,[f1,f2]);
end
	
