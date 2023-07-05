% Tue 21 Sep 14:05:17 CEST 2021
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
%% test a periodogram for hidden periodic frequency components
%% 
%% [p,stat,ratio] = periodogram_test_periodicity_2d(b, nf, bmsk, fmsk, ns)
%%
%% input:
%%	b    (nx * ny): image to test for presence of hidden periodicities,
%%             i.e. periodicities where the frequency is not known a priori
%%      nf   = nfr or [nfx, nfy]
%%	       radius of circular disk (in number of bins) used for smoothing
%%             the periodogram to estimate the spectral density,
%%	       or axes of ellipses for smoothing
%%	       when b is not square a good choice is nfx/nfy ~ Lx/Ly
%%      bmsk : mask in real space selecting parts of the image to include in
%%	       the analysis default is entire image
%%	       the mask can have non-integer values to feather the borders of the mask
%%	fmsk : mask in frequency selecting frequencies to test for periodicity
%%	       default is all frequencies
%%	       note: when b is real, one half plane can always be excluded
%%	       because of symmetry. This slightly increases the significance
%%	ns   : number of samples for the monte-carlo determination of
%%	       the test statistics, mc is only used when parts of the image are masked
%%	       otherwise the analytic test statistic is used 
%%
%% influence of masking the input file:
%% 	      - the root-mean-square energy of the ordinates is proportional
%%	        to the number of unmasked points
%%	      - values in the periodogram are not any more linearly independent
%%	        so that the dof of the filter window is not nf^2
%%
function [pn,stat,out] = periodogram_test_periodicity_2d(b, L, nf, bmsk, fmsk, ns, test_all)
	n = numel(b);
	p_thresh = 0.05;
	if (nargin()< 3 || isempty(L))
		L = n;
	end
	if (nargin()<4)
		bmsk = [];
	end
	if (nargin()<5)
		fmsk = [];
	else
		if (~islogical(fmsk))
			error('fmsk must be logical');
		end
	end	
	if (nargin()<6)
		ns = 100;
	end
	if (nargin<7)
		test_all = true;
	end

	% exclude mean and masked area
	if (~isempty(bmsk))
		b = b-mean(b(bmsk));
                b = b.*bmsk;
	else
		b = b-mean(b(:));
	end
	% compute the 2D periodogram
	Shat  = abs(1/n*fft2(b)).^2;
	% it is not necessary to scale here, as we are testing ratios
	% though the fraction of spectral energy is returned, that is why Shat
	% is normalized to unity over the half-plane
	df    = 1./L;
	Shat  = 2*Shat/(sum(Shat,'all')*df(1)*df(2));

	% estimate the spectral density by smoothing
	% the weights of the smoothing window have to be 1 or zero for the
	% analytic solution of the test statistic,
	% so gaussian smoothing can only be used in combination with the monte-
	% carlo estimating of the statistic, for getting identical results
	% for masking or not masking, we always use a window with weights 0,1 here
	%if (isempty(bmsk))
		[Sbar,nf2] = circfilt2(Shat,nf);
	%else
	%	%note: if this is chosen, the quantiles have to be differently estimated
	%	[Sbar,nf2] = gaussfilt2(Shat,nf);
	%end

	% ratio of periodogram and density
	ratio = Shat./Sbar;

%	% number of tested bins
	if (~isempty(fmsk) && ~isscalar(fmsk))
		nt  = sum(sum(fmsk));
	else
		nt = numel(b);
	end

	% maximum ratio
	if (~isempty(fmsk))
		[ratio_max,mdx] = max(flat(ratio.*fmsk));
	else
		[ratio_max,mdx] = max(flat(ratio));
	end

	if (isempty(bmsk))
		% inclusive
		% Shat/Sbar ~ nf^2*betarnd(a,b)

		% distribution parameters of test statistics
		ap = 1;
		bp = (nf2-1);
		% p-value, if only 1 bin were tested
		p1 = 1-betacdf(ratio_max/nf2,ap,bp);

		% correct for repeated testing the nt-bins selected by the mask
		pn  = 1-(1-p1)^nt;
		if (test_all)
			p1_all  = 1-betacdf(ratio/nf2,ap,bp);
			pn_all  = 1-(1-p1_all).^nt;
		end

		if (0) %nargout()>2)
			np = 100;
			out.pr1 = (1:np)'/(np+1);
			out.qr1 = nf2*betainv(out.pr1,pa,pb);
			out.prn = 1-(1-out.pr1).^nt;
			out.qrn = out.qr1;
		end
	else
		% approximate test statistic by the monte-carlo method
		[pr,qr1,qrn]   = approximate_ratio_distribution(bmsk,nf,ns,fmsk,mdx);
		if (ratio_max > qr1(end))
			% extrapolate
			p1 = 0.5/ns;
		else
			p1  = 1-interp1(qr1,pr,ratio_max,'linear');
		end
		if (ratio_max > qrn(end))
			% extrapolate
			pn = 0.5/ns;
		else
			pn = 1-interp1(qrn,pr,ratio_max,'linear');
		end
		if (test_all)
			p1_all = [];
			% note that while it is defined here for all bins,
			% it is only valid in the range defined by the mask
			pn_all = 1-interp1(qrn,pr,ratio,'linear');
			pn_all(ratio>qrn(end)) = 0.5/ns;
		end
		out.pr1 = pr;
		out.prn = pr;
		out.qr1 = qr1;
		out.qrn = qrn;
	end
	if (test_all)
		fdx = (pn_all<=p_thresh);
		stat.intShat_sig = 0.5*sum(Shat(fdx))*df(1)*df(2);
	end
	[fx,fy]            = fourier_axis_2d(L,size(b));
	[id,jd]            = ind2sub(size(b),mdx);
	stat.max.ratio     = ratio_max;
	stat.max.Shat      = Shat(mdx);
	stat.max.Shat_rel  = Shat(mdx)*df(1)*df(2);
	stat.max.fx        = fx(id);
	stat.max.fy        = fy(jd);
	stat.p1            = p1;
	stat.pn            = pn;

	if (test_all)
		stat.p1_all = p1_all;
		stat.pn_all = pn_all;
	end

	if (nargout()>2)
		out.p = (1:nt)'/(nt+1);
		out.q(:,1) = sort(flat(ratio(fmsk)));

		out.Shat = Shat;
		out.Sbar = Sbar;
		out.ratio = ratio;
	end
end % periodogram_test_periodicity_2d

