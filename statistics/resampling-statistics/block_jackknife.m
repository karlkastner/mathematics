% 2016-04-28 12:53:31.487698856 +0200
% Karl Kastner, Berlin
%
% Jacknfife on blocks for autocorrelated data
function s2 = block_jackknife(fun,x,h)
	n = size(x,1);
	n2 = size(x,2);

	% get estimate for complete set
	Tn  = fun(x);
	k   = floor(n/h);
	Tni = zeros(k,n2);
	% get estimates with deleted blocks
	for idx=1:k
		% mask
		mask = true(n,1);
		mask(1+h*(idx-1):h*idx) = false;
		% get estimate for current subset
		Tni(idx,:) = fun(x(mask,:));
	end
	% pseudo variable
	Tni_ = 1/h*bsxfun(@minus, n*Tn, (n-h)*Tni);
%	s2 = 1/(n*m) * sum( (Tni_ - 1/m*sum(Tni_)).^2 );
	% Carlstein
	% variance
	s2 = h/(k-1) * sum( bsxfun(@minus, Tni_, 1/k*sum(Tni_,1)).^2 ,1 );
	% error variance
	s2 = s2*1/n;
end

