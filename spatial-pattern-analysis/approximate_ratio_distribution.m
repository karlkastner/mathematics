% 2022-09-23 15:54:22.164843714 +0200
%
% input :
% bmsk : region selected for periodicity test (smoothes the periodogram)
% nf   : radius of smoothing window (in bins) for estimating the spectral density
% nsample : number of repetitions to estimate the ratio distribution
%           recommended at
%
% output:
% ratio : ratios for each frequency bin and iteration
%
% note the following complications:                             
%       - problem 1 : ratio locally differs near fr=0, fx,fy=fmax and fx,fy=fmax/2                             
%	- fits of the fisher or beta distribution are highly unstable
%
function [pr, qr1, qrn, ratio, w] = approximate_ratio_distribution(bmsk,nf,ns,fmsk,mdx)
	n     = size(bmsk);

	% white noise (Gaussian field w/o spatial correlation)
	x     = randn(n(1),n(2),ns);

	% periodogram, normalization not necessary, as we determine ratios
	Shat  = abs(fft2(bmsk.*x)).^2;

	% approximate the spectral density by smoothing
	[Sbar, nf2, w,] = circfilt2(Shat, nf);

	ratio = Shat./Sbar;

	% p-values
	pr = (1:ns)'/ns;
	[mid,mjd] = ind2sub(size(bmsk),mdx);

	% associated quantiles of the mth-frequency bin
	qr1 = squeeze(ratio(mid,mjd,:));
	qr1 = sort(qr1);

	% associated quantiels of the maximum of all frequency bins defined by mask msk
	qrn = zeros(ns,1);
	for idx=1:ns
		qrn(idx) = max(flat(ratio(:,:,idx).*fmsk));
	end
	qrn = sort(qrn);
end

