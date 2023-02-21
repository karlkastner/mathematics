% So 19. Jul 12:45:11 CEST 2015
% Karl Kastner, Berlin
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
%% weighted variance of columns, corrected for degrees of freedom (bessel)
%% variance of the weighted sample mean of samples with same mean (but not necessarily same variance)
%% s^2 = sum (w^2(x-sum(wx)^2))
%%
%% s2_mu : error of mean, s2_mu : sd of prediction
%
% function [s2 dof] = wvar(w,x)
%
% f = mu_w = (1/sum w) sum w x
% df/dxi = w_i
% s^2 = sum (df/dxi)^2 s_xi^2 = sum wi^2 s_xi^2
%
function [s2, s2_mu, dof] = wvar(w,x,fullpop,varargin)
	if (isvector(x))
		x = cvec(x);
	end
	if (isvector(w))
		w = cvec(w);
		w = repmat(cvec(w),1,size(x,2));
	end
	wx     = w.*x;

	sw     = sum(w);

	% weighted mean
	mu     = sum(wx)./sw;

	% weighted squared residuals
	w2dx2  = w.*bsxfun(@minus,x,mu).^2;

	s2      = sum(w2dx2)./sw;
%	s2       = dof.*s2;

	if (nargin()<3 || ~fullpop)
		sw2    = sum(w.^2);
		% number of degree of freedoms (normalisation again not necessary)
		dof = sw.^2./sw2;
		s2  = dof./(dof-1).*s2;
	end
end % wvar

