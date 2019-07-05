% Mon 22 Aug 15:48:52 CEST 2016

function test_nlcg_6_order

% polynomial optimisation is problematic if the second derivative of the obj function is zero

c = 2.^(0:1);

meta = delft3d_metadata();
x0   = zeros(length(c),1);
opt = meta.opt;
opt.ls_maxiter = 3;
opt.maxiter = 5;
%opt.lsfun   = @line_search_polynomial;
%opt.lsfun   = @line_search_quadratic;
order = 4;
opt.lsh=[];
opt.lsfun = @(fun,x,f0,g,dir,h,lb,ub,maxiter) line_search_polynomial(fun,x,f0,g,dir,h,lb,ub,maxiter,order)
opt
nlcg(@(x) fun(x,c), x0, opt)

end

function [f g] = fun(x,c)
	p = 6;
	f = c*(cvec(x)-1).^p + 1e-1*c*(cvec(x)-1).^2 + 1;
	%f = c*(cvec(x)-1).^p + c*(x-1).^3 + 1;
	if (nargout() > 1)
		g = grad(@(x) fun(x,c), x);
		g = cvec(g);
	end
end

