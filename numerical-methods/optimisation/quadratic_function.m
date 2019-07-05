% Sa 6. Feb 14:23:42 CET 2016
% Karl Kastner, Berlin
%% evaluate quadratic function in higher dimensions
function f = quadratic_function(H,g,f0,x0)
	f = @(x) f0 + g'*(x-x0) + 0.5*(x-x0)'*H*(x-x0);
end

