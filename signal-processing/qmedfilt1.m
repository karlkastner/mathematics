% Wed 12 Sep 14:24:10 CEST 2018
%% medfilt1, after fitting a quadratic polynomial
function yf_ = qmedfilt1(y,nf)
	if (isvector(y))
		y = cvec(y);
	end
	n = size(y);
	%n  = length(y);
	nf = max(nf,3);
	%o  = floor(nf/2);
	m = floor(nf/2);	

	yf = NaN(n(1),nchoosek(nf,3));

	% TODO, this can be vectorized
	for ddx=1:n(2)
	for ndx=1:n
	k=0;
	% quadratic, all triples
	for idx=-m:m-2
	 for jdx=idx+1:m-1
	  for kdx=jdx+1:m
		if (idx+ndx>0 && kdx+ndx<=n(1))
		% fit
		A = vander_1d([idx;jdx;kdx],2);
		c = A \ y(ndx+[idx;jdx;kdx],ddx);
		% predict % A = [1,0,0];
		k = k+1;
		yf(ndx,k) = c(1);
		end
	  end % kdx
	 end % jdx
	end % idx
	end % for n
	% compute median
	yf_(:,ddx) = nanmedian(yf,2);
	end
end

