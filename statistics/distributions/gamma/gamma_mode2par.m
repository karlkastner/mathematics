% Mon 16 Jan 14:13:39 CET 2023
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
function [a,b] = gamma_mode2par(xm,ym,p0,varargin)
	
	if (nargin()<3||isempty(p0))
		% normal approximation
		[mu,sd] = normpdf_mode2par(xm,ym);
		% now we know that
		[p0(1), p0(2)] = gamma_moment2par(mu,sd);
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
	res       = [xm-xm_; log(ym)-log(ym_+eps)];
end

end

