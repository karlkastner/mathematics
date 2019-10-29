% Tue Jun 19 17:56:07 MSK 2012
% Karl Kästner, Berlin

function P_ = project_circle(P, r0, x0);
	% this kind of projectio is problematic, as it looses accuracy
	% if x0 >> 0 and h_min ~ eps
	P_ = [P(:,1) - x0(1), P(:,2) - x0(2)];
	S = r0./sqrt(sum(P_.^2,2));
	P_(:,1) = P_(:,1).*S + x0(1);
	P_(:,2) = P_(:,2).*S + x0(2);
end % project_circle()

