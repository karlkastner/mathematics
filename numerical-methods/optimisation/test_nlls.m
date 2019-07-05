% Mon  5 Sep 00:16:03 CEST 2016
% Karl Kastner, Berlin
function test_nlls()

m = 5;
n = 10000;
x = rand(n,1);
c0 = [1; 10];
r = 1e-3*randn(size(x));

c = [0; 0];
C = [];
R = [];

opt.reltol = 1e-15;
[c res R g gn A] = nlls(@fun,c,opt);
c
R
gn
R = R/(R(1));
gn = bsxfun(@times,gn,1./(gn(1)));
semilogy(1:length(R),[R gn])
%1/n*A'*A



function [res g] = fun(c)

	res = fun_(c) - (fun_(c0)+r);

%%	f = w'*(p*abs(x-1).^4) + w'*(1-p)*abs(x-1).^2;
%	f = w'*(x-1).^2;
%	f = w'*(p*abs(x-1).^4) + (max(w)+1-w)'*(1-p)*abs(x-1).^2;
%	f = w'*(abs(x-1).^p);
%	f = f + w'*((x-1).^2);
	if (nargout() > 1)
		g = grad(@fun,c,[],'one-sided');
	end
end


function f = fun_(c)
	f = x*c(1) + x*c(1).^2 + c(2)*sin(2*pi*x) + c(2).^2*sin(2*pi*x);
end

end
