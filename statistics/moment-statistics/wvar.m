% So 19. Jul 12:45:11 CEST 2015
% Karl Kastner, Berlin
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

	p      = 1;

	sw     = sum(w);
	swp    = sum(w.^p);

	% weighted mean
	mu     = sum(wx)./sw;

	% weighted squared residuals
	w2dx2  = w.^p.*bsxfun(@minus,x,mu).^2;

	s2      = sum(w2dx2)./swp;
%	s2       = dof.*s2;

	if (nargin()<3 || ~fullpop)
		sw2    = sum(w.^2);
		% number of degree of freedoms (normalisation again not necessary)
		dof = sw.^2./sw2;
		s2  = dof./(dof-1).*s2;
	end
end % wvar

