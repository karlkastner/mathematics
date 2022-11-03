% Fri 10 Dec 11:22:33 CET 2021
% Karl Kästner, Berlin
%
%% filter with a rectangular window along both dimensions
% function z = meanfilt2(z,nf)
function xf = meanfilt2(x,nf)
	if (length(nf)<2)
		nf = [nf(1),nf(1)];
	end
	xf = zeros(size(x));
	for idx=1:size(x,3)
		xi = padd2(x(:,:,idx),nf);
		xi = meanfilt1(xi,nf(1));
		xi = meanfilt1(xi.',nf(2)).';
		xf(:,:,idx) = xi(nf(1)+1:end-nf(1),nf(2)+1:end-nf(2));
	end
end
