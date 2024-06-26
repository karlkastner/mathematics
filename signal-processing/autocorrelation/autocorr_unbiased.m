% Wed 26 Jun 09:40:23 CEST 2024
% Karl Kastner, Berlin
%% autocorrelation, truncated at the end and only computed within the interior
%% of the domain, this estimated is unbiased
function a = autocorr_unbiased(y)
	if (isvector(y))
		y = cvec(y);
	end
	siz = size(y);
	a = zeros(siz);
	y = y-mean(y);
	for jdx=1:siz(1)
		a(jdx,:) = mean(y(1:end-jdx+1,:).*y(jdx:end,:));
	end
	% normalize with variance
	a = a./a(1,:);
end

