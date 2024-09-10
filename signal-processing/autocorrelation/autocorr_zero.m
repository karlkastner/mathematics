% 2024-06-26 09:42:13.552258584 +0200
% Karl Kastner, Berlin
%% autocorrelation, values beyond domain end are padded with 0
%% this yields the same result as the autocorr function in the econometric
%% toolbox of matlab 
%% autocorr_zeros(y) = autocorr(y)
function a = autocorr_zero(y)
	if (isvector(y))
		y = cvec(y);
	end
	n = size(y,1);
	y = y-mean(y);
	y = [y; zeros(size(y))];
	a = ifft(abs(fft(y)).^2);
	a = a./a(1,:);
	% clip padded zeros
	a = a(1:n,:);
end

