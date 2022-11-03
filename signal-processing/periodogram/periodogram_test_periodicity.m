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
function [p,ratio,maxShat,mdx,fdx,S] = periodogram_test_periodicity( ...
					       fx,Shat,nf,fmin,fmax,S,mode)
	if (nargin()<4 || isempty(fmin))
		fmin = 0;
	end
	if (nargin()<5 || isempty(fmax))
		fmax = inf;
	end
	if (nargin()<7 || isempty(mode))
		mode = 'inclusive';
	end

	% select frequency range
	fdx = find(fx > fmin & fx < fmax);
	% number of tested bins
	m   = length(fdx);

	% estimate the spectral density
	if (nargin()>6 && ~isempty(S))
		% use provided exact spectral density
		mode = 'exact';
	else
		if (0==mod(nf,2))
			% make odd
			nf = nf+1;	
		end
		% smoothing window
		w = ones(nf,1);
	    	switch (mode)
			case {'exclusive'}
				w((nf+1)/2) = 0; 
			case {'inclusive'}
				% nothing to do
	   	 otherwise
			error('here')
	    	end
		% normalize
		w  = w/sum(w);
		% smooth the periodogram
		S = conv(Shat,w,'same');
	end % else of ~isempty(S)

	% ratio = Shat./Sf;
	[ratio,mdx] = max(Shat(fdx)./S(fdx));
	mdx         = fdx(mdx);
	maxShat     = Shat(mdx);
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
		% Shat/Sf ~ nf*betarnd(a,b)
		a = 1;
		b = (nf-1);
		p1 = 1-betacdf(ratio/nf,a,b);
		% correct for repeated testing
		p  = 1-(1-p1)^m;
	end % switch mode
end % periodogram_test_periodicity

