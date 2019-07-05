% Mon Jan 12 10:38:15 CET 2015
% Karl Kastner, Berlin
%
%% mean angle
% TODO this is not good, use mean of sin and cos
function	alpha = meanangle(alpha)
	fdx = alpha > pi;
	alpha_      = alpha;
	alpha_(fdx) = 2*pi-alpha_(fdx);

	if (var(alpha) < var(alpha_))
		alpha = mean(alpha);
	else	
		alpha = mean(alpha_); %-pi)+pi;
	end
end

