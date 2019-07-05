% Mi 8. Apr 12:02:49 CEST 2015
% Karl Kastner, Berlin
%
%% pooled mean and standard deviation of several groups of different size, mean and standard deviation
function [mu, s, n] = pool(mu,s,n,serrflag)
%	mu_ = mu(1,:);
%	s_  = s(1,:);
%	n_  = n(1,:);
%	% iterative computation
%	for idx=2:size(mu,1)
%		[mu_ s_ n_] = pool_(mu_,s_,n_,mu(idx,:),s(idx,:),n(idx,:));
%	end
	% s is standard error not standard deviation
	if (nargin() > 3 && serrflag)
		s = sqrt(s.*s.*n);
	end

	% computation by recursion, this is numerically more stable than
	% iteration, but takes twice as many steps then iteration
	nn = size(mu,1);
	if (nn > 1)
		% devide
		nn = floor(nn/2);
		% recurse
		[mu1 s1 n1] = pool(mu(1:nn,:),s(1:nn,:),n(1:nn,:));
		[mu2 s2 n2] = pool(mu(nn+1:end,:),s(nn+1:end,:),n(nn+1:end,:));
		% merge
		[mu s n] = pool2(mu1,s1,n1,mu2,s2,n2);
	end
	if (nargin() > 3 && serrflag)
		% convert standard deviation to standard error
		s = sqrt(s.*s./n);
	end
end
% for two groups
% TODO this is not yet bessel corrected for finite samples sizes
function [mu s n] = pool2(mu1,s1,n1,mu2,s2,n2)
	n  = n1+n2;
	mu = (n1.*mu1+n2.*mu2)./(n1+n2);
	ss = (n1.*s1.*s1 + n1.*n2./(n1+n2).*(mu1-mu2).^2 + n2.*s2.*s2)./(n1+n2);
	s  = sqrt(ss);
end

