% 2016-02-27 15:46:30.584686222 +0100
% Karl Kastner, Berlin
%
%% low-pass filter of data
% function f=lowpass2(n,fs,f1,winstr)
function f=lowpass2(n,fs,f1,winstr)
	if (nargin() < 4)
		winstr = 'kaiserwin';
	end

	% normalise sample frequencies
	f1  = f1/fs;

	% design the filter
	f   = fir1(n-1,[f1]);

	% apply a window to the filter kernel
%	win = hanwin(0:n,0.5*(n),0.5*(n));
	if (~isempty(winstr))
		win = feval(winstr,0:n-1); %,0.5*n,0.5*n);
		f   = f.*win;
	end
	% normalise the filter
	f   = f/sum(f);
end

