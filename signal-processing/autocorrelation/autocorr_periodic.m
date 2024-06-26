% 2024-06-26 09:42:13.552258584 +0200
% Karl Kastner, Berlin
%% autocorrelation, values beyond domain end periodically extended
function a = autocorr_periodic(y)
	y = y-mean(y);
	a = ifft(abs(fft(y)).^2);
	a = a./a(1,:);
end
