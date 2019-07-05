% 2016-03-07 16:11:47.075917226 +0100
%
%% proxy standard deviation determined by quantiles
function sd = qstdq(q,p)
	if (isvector(q))
		q = cvec(q);
	end
	scale = 1/(2*norminv(1-p));
	sd = scale*(q(2,:)-q(1,:));
end

