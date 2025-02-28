% Fri 31 Jan 11:52:20 CET 2025
function rho = corr1d(z)
	% zij = 1/2*(zl + zr)*rho + eps
	% ((1 - r/2)I - r/2*D) z = eps
	rho = zeros(size(z,2),1);
	z = z-mean(z);
	for idx=1:size(z,2)
		rhs      = (up(z(:,idx),1) + down(z(:,idx),1));
		%rho(idx,1) = (z(:,idx) \ rhs(:))/2;
		rho(idx,2) = (z(:,idx)'*rhs)/(z(:,idx)'*z(:,idx))/2;
	end
end

