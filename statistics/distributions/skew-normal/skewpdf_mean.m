% 2024-05-15 10:07:38.194747129 +0200
% Karl Kastner, Berlin
function mu = skewpdf_mean(a)
	[sk, sd, mu] = skewpdf_central_moments(a);
end

