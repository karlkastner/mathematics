% 2016-03-10 15:04:44.019300629 +0100
%
%% expected norm of x.^n, when values x in x are iid normal with mu and sigma
%
function m = normmoment(mu,sigma,order)
	% TODO use generalised hermite polynomials
	switch (order)
	case {0}
		m = 1;
	case {1}
		m = mu;
	case {2}
		m = mu.^2+sigma.^2;
	case {3}
		m = mu.*(mu.^2+3*sigma.^2);
	case {4}
		m = mu.^4+6*mu.^2.*sigma.^2+3*sigma.^4;
	end
end

