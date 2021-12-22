% Tue 21 Sep 14:05:17 CEST 2021
% mode : exact, inclusive, exclusive
% function [p,ratio,maxShat,mdx,fdx,S] = periodogram_test(fx,Shat,fmin,fmax,mode,varargin)
function [p,ratio,maxShat,mdx,fdx,S] = periodogram_test(fx,Shat,fmin,fmax,mode,varargin)
	if (nargin()<3)
		fmin = 0;
	end
	if (nargin()<4)
		fmax = inf;
	end
	if (nargin()<5)
		mode = 'inclusive';
	end

	fdx = find(fx > fmin & fx < fmax);
	m   = length(fdx);

	switch (mode)
	case {'exact'}
		S = varargin{1};
	case {'exclusive'}
		nf = varargin{1};
		if (0==mod(nf,2))
			% make odd
			nf = nf+1;	
		end
		w = ones(nf,1);
		w((nf+1)/2) = 0; 
		w  = w/sum(w);
		% smooth the periodogram
		S = conv(Shat,w,'same');
	case {'inclusive'}
		nf = varargin{1};
		if (0==mod(nf,2))
			% make odd
			nf = nf+1;	
		end
		w = ones(nf,1);
		w = w/sum(w);
		% smooth the periodogram
		S = conv(Shat,w,'same');
	otherwise
		error('here')
	end

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
	end
end

