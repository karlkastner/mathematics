% So 19. Jul 12:45:11 CEST 2015
% Karl Kastner, Berlin
%
%% weighed standard deviation
%
function [sd] = wstd(w,x)
	sd = sqrt(wvar(w,x));
%	if (isvector(x))
%		x = cvec(x);
%		w = cvec(w);
%	end
%		wx  = w.*x;
%		fdx = isfinite(wx);
%		mu  = zeros(1,size(x,2));
%		sw  = zeros(1,size(x,2));
%		n2  = size(x,2);
%		for idx=1:n2
%			sw(idx) = sum(w(fdx(:,idx),idx));
%			mu(idx) = sum(wx(fdx(:,idx),idx))/sw(idx);
%		end
%		s2 = zeros(1,size(x,2));
%		wdx2  = w.*bsxfun(@minus,x,mu).^2;
%		for idx=1:n2
%			% number of degree of freedoms
%			ni      = sw(idx).^2./sum(w(fdx(:,idx),idx).^2);
%			% this is finite population corrected
%			s2(idx) = (ni/(sw(idx)*(ni-1))).*sum(wdx2(fdx(:,idx),idx));
%		end
%		sd = sqrt(s2);
end % wstd

