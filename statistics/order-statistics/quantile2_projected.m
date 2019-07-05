% Mo 29. Feb 16:37:25 CET 2016
%
%% quantile in two dimensions
%
function [Q q a] = quantile2_projected(x,p)
%	p = (1+p)/2;
	me = median(x,2);
	n = size(x,2);
	% direction
	d = bsxfun(@minus,x,me);
	no = 1./sqrt(sum(d.^2));
	dir = bsxfun(@times,d,no);
	% sort by angle
	a = atan2(dir(1,:),dir(2,:));
	[a sdx] = sort(a);
	% x does not have to be sorted
	x = x(:,sdx);
	dir = dir(:,sdx);
	d   = d(:,sdx);
	q = zeros(1,n);
	for idx=1:n
		% inner product with direction
%		d_ = dir(:,idx)'*d;
%		d_ = d_.*d_;
%		w = quantile(d_,0.84);
%		d_ = (dir(:,idx)*dir(:,idx)')*d;
%		Q(:,idx) = quantile(d_',p)';
		

%		fun = @(w) sum(dir(:,idx)'*d < w) - n*p;
%		w = 0;
%		while (fun(w) < 0)
%			w = w+0.1;
%		end

%		for jdx=0:50
%			fun(jdx/10)
%		end
%		w = lsqnonlin(@(w) sum(dir(:,idx)'*d < w) - n*(2*p-1), 1);
%		Q(:,idx) = w*dir(:,idx);

		w = dir(:,idx)'*d;
		w = w.*no;
%		w = bsxfun(@times,dir(:,idx),d);
%		w = sign(w).*(abs(w));
%		w = sum(w);
		%w = d(:,idx)'*dir;
		%w = dir(:,idx)'*dir;
		q(idx)   = quantile(w,p);
%		w = sign(q(idx)).*q(idx)^0.5;
		Q(:,idx) = me + q(idx)*dir(:,idx);
		%Q(:,idx) = me + q(idx)*d(:,idx);
		%P = eye(2) - dir(:,idx)*dir(:,idx)';
		%P = dir(:,idx)*dir(:,idx)';
		%xp = P*d;
		%Q(:,idx) = me+quantile(xp',p)';
	end
	%Q =bsxfun(@plus,me,bsxfun(@times,q,dir));
end
