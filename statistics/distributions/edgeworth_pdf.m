% Do 3. Mär 21:06:44 CET 2016
% Karl Kastner, Berlin
%%
%% probability density of and unknown distribution
%% with mean mu, standard deviation sigma, and third and fourth cumulants
%% c.f. Rao 2010
function f = edgeworth_pdf(x,mu,s,sk,ku,n)
	% normalise x
	x = (x-mu)/s;
	H = @(n,x) hermiteH(n,x);
	f = normpdf(x).*(  1 ... 
			 + 1/6*sk*H(3,x) ...
			 + 1/24*(ku-3)*H(4,x) ...
			 + 1/72*sk^2*H(6,x));
end

