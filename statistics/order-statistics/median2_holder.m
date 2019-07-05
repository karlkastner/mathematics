% 2016-03-01 09:16:31.349392965 +0100
%
% two-dimensional median
% c.f. bedall and zimmerman
function me = median2_holder(X)
	me = rand(2,1);
	me = lsqnonlin(@(q) objective(X,q),me);
end
function f = objective(X,q)
	 f = sum(bsxfun(@minus,X,q).^2,2);
end


