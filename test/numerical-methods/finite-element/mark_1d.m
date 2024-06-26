% Mon Apr 30 12:12:57 MSK 2012
% Karl Kästner, Berlin

% P : set of points, including boundary
% V : values at the points
% T : elements (line segments)
function [M err_est v_err] = mark_1d(P, T, V, N, h_side, C)

	lp = size(P,1);
	lt = size(T,1);

	% allocate memory
	v_err = zeros(lt,1);

	% calculate first derivative inside the elements
	dV = (V(T(:,1)) - V(T(:,2))) ./ (P(T(:,1)) - P(T(:,2)));

%	% calculate element centres
%	C = 0.5*(P(T(:,1)) + P(T(:,2)));

	% calculate estimated norm of the second derivative per element
	for idx=1:lt
		s = 0;
		for jdx=1:2
		 if (N(idx,jdx) > 0)
			a = idx;
			b = N(idx,jdx);
			s = max(s, abs((dV(a) - dV(b))./(C(a) - C(b))));
		 end
		end
		% d^2 u/dx * len
		v_err(idx,1) = s * h_side(idx).^2; %abs( P(T(idx,1)) - P(T(idx,2))).^2;
	end

	% get error estimate
	err_est = max(v_err);

	% mark elements for refinement	
	M = find(v_err > 0.5^2*err_est);
end % mark_1d()

