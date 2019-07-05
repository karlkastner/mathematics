% Do 18. Feb 14:49:13 CET 2016
% Karl Kastner, Berlin
%
%% two dimensional oja median
%% note: the multivariate median is not unique
%%
%% oja 1983, for extension to multivariate function, see chaudhri
function me = oja_median2(x)
	x0 = mean(x,2);
	me = lsqnonlin(@(x0) delta(x,x0), x0);
end

function s = delta(x,x0)
	% compute distance for each pair
	n = size(x,2);
	A = ones(3);
	A(2:3,end) = x0;
	s = 0;
	for idx=1:n-1
		A(1,2:3) = x(:,idx);
		for jdx=idx+1:n
			A(2,2:3) = x(:,jdx);
			s = s+0.5*abs(det(A));
		end
	end
end


