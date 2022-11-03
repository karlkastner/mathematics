% 2016-07-09 16:52:14.382477863 +0200
% Karl Kastner, Berlin
%
%% design low pass filter with cutoff-frequency f1
%
% function y  = lowpass2(x,n,fs,f1,nanflag,exflag);
% TODO treatment of boundary
function y  = lowpass2(x,n,fs,f1,nanflag,exflag);
	if (nargin() < 5)
		nanflag = 0;
	end
	if (nargin() < 6)
		exflag = 0;
	end

	win = lowpassfilter2(n,fs,f1);
	if (nanflag)
		fdx = isnan(x);
		x(fdx) = 0;
		o = ones(size(x));
		o(fdx) = 0;
		s = conv_centered(o,win);
		s(s<0.5) = NaN;
	end

	if (exflag)
		length(x)
		x_ = x;
		x   = [x(1)*ones(n,1);
                       cvec(x);
		       x(end)*ones(n,1)];
		y   = conv_centered(x,win);
		y   = y(n+1:end-n);
	else
		x   = [NaN;cvec(x);NaN];
		y   = conv_centered(x,win); %conv(x,win,'same');
		y   = y(2:end-1);
	end
	if (nanflag)
		y = y./s;
	end
end

