% Tue 21 Sep 14:05:17 CEST 2021
% Karl Kästner, Berlin
%
%% test a periodogram for hidden periodic frequency components
%% 
%% [p,stat,ratio] = periodogram_test_periodicity_2d(b, L, nf, frmin, frmax)

%% input:
%%	fx : frequengcies
%%	b    : image to test for presence of hidden periodicities,
%%             i.e. periodicities where the frequency is not known a priori
%%      nf   : radius of circular disk (in number of bins) used for smoothing
%%             the periodogram to estimate the spectral density
%%      bmsk : mask determining parts of the image to include in the analysis
%%             default is entire image
%%	fmin, fmax : (radial) frequency range limits to test (fmask)
%%
%%	mode : automatically set to "exact", when S given
%%	       inclusive : estimate density by smoothing including the central bin
%%	       exclusive : estimate density by smoothing  excluding the central bin
%%	note: inclusive and exclusive lead to different distribution
%%	      but identical p-values
%%
%% influence of masking the input file:
%% 	      - the root-mean-square energy of the ordinates is proportional
%%	        to the number of unmasked points
%%	      - values in the periodogram are not any more linearly independent
%%	        so that the dof of the filter window is not nf^2
%%
%% TODO make frmin, frmax an fmask
function [p,stat,ratio,qq] = periodogram_test_periodicity_2d(b, L, nf, bmsk, fmsk, ns)
	mode = 'inclusive';
	if (isempty(L))
		L = [1,1];
	end
%	if (length(nf)<2)
%		nf = [nf,nf];
%	end
	if (nargin()<4)
		bmsk = [];
	end
	if (nargin()<5)
		fmsk = [];
	end
%	if (nargin()<5||isempty(frmin))
%		frmin = 0;
%	end
%	if (nargin()<6||isempty(frmax))
%		frmax = inf;
%	end
	n = size(b);	
	[fx, fy, fr, ft] = fourier_axis_2d(L,n);

	n = numel(b);

	% exclude mean and masked area
	if (~isempty(bmsk))
		b = b-mean(b(bmsk));
		b(~bmsk) = 0;
	else
		b = b-mean(b(:));
	end
	% compute the 2D periodogram
	% it is not necessary to scale here, as we are testing ratios
	Shat  = abs(1/n*fft2(b)).^2;

	% estimate the spectral density by smoothing
	% S     = ifftshift(meanfilt2(fftshift(Shat), nf));
	if (isempty(bmsk))
		[Sbar,nf2] = circfilt2(Shat,nf);
	else
		[Sbar,nf2] = gaussfilt2(Shat,nf);
	end

	% ratio of periodogram and density
	ratio = Shat./Sbar;

%	% select frequency range
%	fmsk = (fx >= 0) & (fr >= frmin) & (fr <= frmax);

%	% number of tested bins
	nt  = sum(sum(fmsk));

	% maximum ratio
	[ratio_max,mdx] = max(flat(ratio.*fmsk));

	if (isempty(bmsk))
		% inclusive
		% Shat/Sbar ~ nf^2*betarnd(a,b)
		a = 1;
		b = (nf2-1);
		p1 = 1-betacdf(ratio_max/nf2,a,b);
		% correct for repeated testing
		p  = 1-(1-p1)^nt;
	else
		[pr,qr1,qrn]   = approximate_ratio_distribution(bmsk,nf,ns,fmsk,mdx);
		if (ratio_max > qr1(end))
			% extrapolate
			p1 = 0.5/ns;
		else
			p1  = 1-interp1(qr1,pr,ratio_max,'linear');
		end
		if (ratio_max > qrn(end))
			% extrapolate
			p = 0.5/ns;
		else
			p = 1-interp1(qrn,pr,ratio_max,'linear');
		end
	end

	stat.ratio_max = ratio_max;
	stat.fr_max    = fr(mdx);
	% stat.Shat_max  = Shat(mdx); this is not passed, as unscaled
	stat.p1 = p1;

	if (nargout()>3)
		qq.p = (1:nt)'/(nt+1);
		qq.q(:,1) = sort(flat(ratio(fmsk)));
		%qq.q(:,2) = nf2*betainv(qq.p,a,b);
	end
end % periodogram_test_periodicity_2d

