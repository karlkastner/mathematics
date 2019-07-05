% 2013-08-17 21:41:42.384959347 +0200
%
%% nearest neighbour interpolation in two dimensions
% Ip : nearest neighbour
function [Vp Ip] = interp2_man(X,Y, V,Xp,Yp)
	XX = X(:)*ones(1,length(Y));
	YY = ones(length(X),1)*Y(:)';
	for idx=1:prod(size(Xp))
		[dist2 Ip(idx)] = min( (XX(:)-Xp(idx)).^2 + (YY(:)-Yp(idx)).^2 );
		Vp(idx,1) = V(Ip(idx));
	end
end % interp2_man

