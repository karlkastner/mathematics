% 2015-10-12 18:52:33.639617864 +0200 binormrnd.m
% Karl Kastner, Berlin
% 
%% generate two correlated normally distributed vectors
% TODO call this as subfunction of randar1
function [X, Y] = binormrnd(mu1,mu2,s1,s2,rho,n)
	s1 = s1.^0.5;
	s2 = s2.^0.5;
	if (abs(rho) > 1)
		error('')
	elseif (rho==1)
		C = [s1,+s2; 0 0];
	elseif (rho == -1)
		C = [s1,-s2; 0 0];
	else
		C=[s1^2, rho*s1*s2;rho*s1*s2,s2.^2];
		C=chol(C);
	end
	X = randn(n,2)*C + ones(n,1)*[mu1 mu2];
	if (nargout > 1)
	Y = X(:,2);
	X = X(:,1);
	end
end

