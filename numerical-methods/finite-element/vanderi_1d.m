% Fri  6 Jan 15:46:15 CET 2017
% Karl Kästner, Berlin
%
% vander monde matrix of the dth-derivative
%function V = vanderd_1d(x,d,order)
% TODO change definition to make d last argument
function V = vanderi_1d(x,d,order)
	for idx=1:order+1
		V(:,idx) = factorial(idx-1)./factorial(idx+d-1)*x.^(idx-1+d);
	end

%	if (0 == d)
%		% no derivative
%		V = vander_1d(x,order);
%	else
%		V      = zeros(size(x,1),order+1,class(x));
%		if (d<=order)
%		V(:,d+1) = factorial(d);
%	%	for idx=d+1:order+1
	%		V(:,idx) = V(:,idx-1).*x;
	%	end
	%	for idx=d:order+1
	%		V(:,idx) = factorial(idx)/factorial(idx-d)*V(:,idx);
	%	end
%		for idx=d+2:order+1
%			V(:,idx) = (idx-1)/(idx-d-1)*V(:,idx-1).*(x);
%		end
%		end
%	end
end % vanderi_1d

