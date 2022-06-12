% Tue 21 Sep 14:05:17 CEST 2021
% Karl Kästner, Berlin
%
%% test a periodogram for hidden periodic frequency components
%% 
%% function [p,ratio,maxShat,mdx,fdx,S] = periodogram_test_periodicity(fx,Shat,nf,fmin,fmax,S,mode)

%% input:
%%	fx : frequengcies
%%	Shat : corresponding periodogram values
%%      nf   : number of bins to test for periodicity, ignored when S is given
%%	fmin, fmax : frequency range limits to test
%%	S    : exact (a priori known theoretical spectral density, must not be estimated from the periodogram)
%%	mode : automatically set to "exact", when S given
%%	       inclusive : estimate density by smoothing including the central bin
%%	       exclusive : estimate density by smoothing  excluding the central bin
%%	note: inclusive and exclusive lead to different distribution
%%	      but identical p-values
function [p,stat,ratio] = periodogram_test_periodicity_2d( b, L, nf, frmin, frmax)
	mode = 'inclusive';
	if (isempty(L))
		L = [1,1];
	end
	if (nargin()<4)
		frmin = 0;
	end
	if (nargin()<5)
		frmax = inf;
	end
	n = size(b);
	[fx, fy, fr, ft] = fourier_axis_2d(L,n);
	
	% compute the 2D periodogram
	% it is not necessary to scale here, as we are testing ratios
	n = numel(b);
	Shat  = abs(1/n*fft2(b-mean(b(:)))).^2;
	%Shat  = fftshift(Shat);
	% estimate the spectral density by smoothing
	S     = ifftshift(meanfilt2(fftshift(Shat), nf));
%	Shat  = Shat(nf:end-nf+1,nf:end-nf+1);
%	S     = S(nf:end-nf+1,nf:end-nf+1);
	ratio = Shat./S;

%	% select frequency range
%	TODO, select only 1/2 of the plane
%	fdx = find(fx > fmin & fx < fmax);
%	% number of tested bins
	fdx = (fx >= 0) & (fr >= frmin) & (fr <= frmax);
	nt  = sum(sum(fdx));

	[ratio_max,mdx] = max(flat(ratio.*fdx));
	%mdx            = fdx(mdx);
	%q           = maxShat/S(mdx);
	switch (mode)
	case {'exact'}
		p = periodogram_p_value(maxShat,S(mdx),m);
	case {'exclusive'}
		% Shat/Sf ~ frnd(d1,d2)
		d1 = 2;
		d2 = 2*(nf-1);
		p1 = 1-fcdf(ratio,d1,d2);
		% correct for repeated testing
		p  = 1-(1-p1)^m;
	case {'inclusive'}
		% Shat/S ~ nf^2*betarnd(a,b)
		a = 1;
		b = (nf^2-1);
		p1 = 1-betacdf(ratio_max/nf^2,a,b);
		% correct for repeated testing
		p  = 1-(1-p1)^nt;
	end % switch mode
	stat.ratio_max = ratio_max;
	stat.fr_max    = fr(mdx);
	% stat.Shat_max  = Shat(mdx); this is not passed, as unscaled
	stat.p1 = p1;

end % periodogram_test_periodicity

