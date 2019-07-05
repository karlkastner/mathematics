% 2016-02-28 00:24:11.066546441 +0100
%% design finite impulse response filter by the least squares method
% M   : length of filter in taps
% fs  : sample frequency
% fun : function returning desired frequency response
function f = firls_man(M,fun,fs)
	% sample interval
	L = round(1/fs);
	% compute colums of the regression matrix
	k     = (0:L)';
	omega = 2*pi*k/L;
	for idx=1:M-1
		n = idx;
		F(:,idx) = 2*cos(omega*(M-n));
	end
	F(:,M) = 1;
	% right hand side
	% where does this sqrt come from?
	a = fun(omega);
	% solve least squares problem
	f = F\a;
	% complete
	% TODO where does the 2 come from?
	f = 2*[f;f(end-1:-1:1)];
end

