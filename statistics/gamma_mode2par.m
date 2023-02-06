% Mon 16 Jan 14:13:39 CET 2023
function [a,b] = gamma_mode2par(xm,ym,p0,varargin)
	
	if (nargin()<3||isempty(p0))
		p0 = [1.01,0.1];
	end
	opt = optimset('MaxFunEvals',1e4,'MaxIter',1e3,'Display','notify');
	l = lsqnonlin(@resfun,p0,[],[],opt);
%	if (flag ~= 1)
%		warning(sprintf('lsqnonlin did terminate without converging %d\n',flag));
%	end
	a = l(1);
	b = l(2);

	% it is not clear to mear, why the univariate optimization performs worse
	% and often does not converge, maybe bc of the condition that the maximum 
	%b0 = 1./ym;
	%b0 = 1.1;
	%[b(2)] = fzero(@(b) gampdf(xm,1+ym.*b,b)-ym,b0);
	%[a(2)] = fzero(@(a) gampdf(xm,a,ym./(a-1))-ym,b0);
	%[b] = fzero(@(b) log(gampdf(xm,1+ym./b,b))-log(ym),b0);
	%a(2)= lsqnonlin(@(b) log(gampdf(xm,a,ym./(a-1)))-log(ym),b0);

function res = resfun(l)
	[xm_,ym_] = gamma_mode(l(1),l(2),varargin{:});
	res       = [xm-xm_; log(ym)-log(ym_)];
end

end

