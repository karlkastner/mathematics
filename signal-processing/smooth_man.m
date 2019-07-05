% 2015-01-21 17:01:08.362756339 +0100
function	Y = smooth(Y,order,cl)
	X = (1:size(Y,1))';
	X = X-mean(X); % length(Y)/2
	X = X/max(X);
%	A = vander_1d(X,order);
	A = fmat([],order,X);
	A = orth(A);
	for idx=1:size(Y,2)
		fdx = find(isfinite(Y(:,idx)));
		c = A(fdx,:) \ Y(fdx,idx);
		if (nargin() > 2)
			res  = Y(fdx,idx) - A(fdx,:)*c;
			res2 = res.*res;
			s2   = quantile(res2,[0.68]);
			gdx  = (res2 < cl^2*s2);
			c = A(fdx(gdx),:) \ Y(fdx(gdx),idx);
		end
		Y(:,idx) = A*c;
	end
end

