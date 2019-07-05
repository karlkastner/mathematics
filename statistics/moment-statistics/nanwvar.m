% So 19. Jul 12:45:11 CEST 2015
% Karl Kastner, Berlin
%
%% weighted variance of columns, corrected for degrees of freedom (bessel)
%%
%% s^2 = sum(w*(x-sum(wx)/sum(w))^2)/sum(w)
%
% function [s2 dof] = wvar(w,x)
%
function [s2 dof] = nanwvar(w,x)
	if (isvector(x))
		x = cvec(x);
		w = cvec(w);
	end
	if (isvector(w))
		w = repmat(cvec(w),1,size(x,2));
	end

	wx      = w.*x;
	fdx     = isnan(wx);
	w(fdx)  = 0;
	wx(fdx) = 0;

	sw    = sum(w);
	sw2   = sum(w.^2);
	mu    = sum(wx)./sw;
	wdx2  = w.*bsxfun(@minus,x,mu).^2;
	wdx2(isnan(fdx)) = 0;

	% number of degree of freedoms (normalisation again not necessary)
	dof = sw.^2./sw2;
	s2  = dof./(dof-1).*(sum(wdx2)./sw);
end % wvar

