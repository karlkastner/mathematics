% Fri Dec 28 22:40:45 MSK 2012
% Karl Kästner, Berlin

function p = cg_coef_to_poly(a,b)
	r = zeros(1,length(a)+1);
	p = zeros(1,length(a)+1);
	r(end) = 1;
	p(end) = 1;

	% todo optimise indices
	for idx=1:length(a)
		r(1:end-1) = r(1:end-1) - a(idx)*p(2:end);
		p = r + b(idx)*p;
	end
end

