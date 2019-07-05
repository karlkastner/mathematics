% 2016-07-09 12:47:47.727971477 +0200
%
%% mean angle
% TODO this is inefficient, simply use sin/cos
% note: if angles are exactly 180deg appart, than the average angle is undefined
function mu = average_angle(a)
	mu = mean(a);
	for idx=1:size(a,2)
		mu(idx) = lsqnonlin(@(mu) objective(a-mu),mu(idx));
	end
end

function d	= objective(d)
		flag    = d < -pi;
		d(flag) = 2*pi + d(flag);
		flag    = d > pi;
		d(flag) = 2*pi-d(flag);
end

