% 2023-01-17 14:51:18.484870247 +0100
% function h = lognpdf_entropy(lmu,lsd)
%% c.f. queiroz 2016
function h = lognpdf_entropy(lmu,lsd)
	if (issym(lmu) || issym(lsd))
		pi_ = sym(pi);
	else
		pi_ = pi;
	end
	h = log2(sqrt(2*pi_)*lsd*exp(lmu+0.5));
end

