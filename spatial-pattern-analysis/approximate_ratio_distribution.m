% 2022-09-23 15:54:22.164843714 +0200
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% input :
%% bmsk : region selected for periodicity test (smoothes the periodogram)
%% nf   : radius of smoothing window (in bins) for estimating the spectral density
%% nsample : number of repetitions to estimate the ratio distribution
%%           recommended at
%%
%% output:
%%	pr    : probabilities for quantiles
%%	qr1   : quantiles of the distribution for bin m
%%	qrn   : quantiles of the distribution for the maximum of bins selected by fmsk 
%% 	ratio : ratios for each frequency bin and iteration (only for last block, for testing)
%% intput: 
%%	bmsk : mask region pattern/interest in the real domain
%%	nf   : smoothing window radius in the frequency domain for density estimation
%%	ns   : number of samples for the monte-carlo simulation
%%	fmsk : mask frequencies of interest 
%%	mdx  : selection of an a-priori known frequency bin
%%
%% note the following complications:                        
%%       - problem 1 : ratio locally differs near fr=0, fx,fy=fmax and fx,fy=fmax/2                             
%%	- fits of the fisher or beta distribution are highly unstable
%%
function [pr, qr1, qrn, ratio, w] = approximate_ratio_distribution(bmsk,nf,ns,fmsk,mdx)
	if (isempty(fmsk))
		fmsk = 1;
	end
	sb     = size(bmsk);

	% allocate memory
	qr1 = zeros(ns,1);
	qrn = zeros(ns,1);

	% monte-carlo simulation of quantiles,
	% avoids running out of memory
	% and is considerably faster than computing it blockwise
	% or event for all samples simultaneously
	for idx=1:ns
		approximate_ratio_distribution_(idx,idx);
	end	

	% probability	
	pr = ((0:ns+1)')/(ns+1);

	% quantiles
	qr1 = sort(qr1);
	qrn = sort(qrn);

	% extrapolation to 0 and 1
	qr1 = [0; qr1; 2*qr1(end)-qr1(end-1)];
	qrn = [0; qrn; 2*qrn(end)-qrn(end-1)];

function approximate_ratio_distribution_(k1,k2)
	% white noise (Gaussian field w/o spatial correlation)
	x     = randn(sb(1),sb(2),k2-k1+1);

	% periodogram, normalization not necessary, as ratios are tested
	% note that we do not have to subtract the mean here, as the expected mean of x is zero
	Shat  = abs(fft2(bmsk.*x)).^2;

	% approximate the spectral density by smoothing
	[Sbar, nf2, w,] = circfilt2(Shat, nf);

	ratio = Shat./Sbar;

	% p-values
	[mid,mjd] = ind2sub(sb,mdx);

	% associated quantiles of the mth-frequency bin only
	qr1(k1:k2) = squeeze(ratio(mid,mjd,:));

	% associated quantiles of the maximum ratio of all frequency bins selected
	% by the frequency mask
	for idx=k1:k2
		qrn(idx) = max(ratio(:,:,idx-k1+1).*fmsk,[],'all');
	end
end % approximate_ratio_distribution_

end % approximate_ratio_distribution

