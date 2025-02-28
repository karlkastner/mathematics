% Fri 31 Jan 12:04:54 CET 2025
% spatial correlation in two dimensions
% assumes circular bc
function rho = corr2d(z)
	% zij = 1/4*(zl + zr + zu + zd)*rho + eps
	% ((1 - r/4)I - r/4*D) z = eps
	rho = zeros(size(z,3),1);
	for idx=1:size(z,3)
		rhs      = (   up(z(:,:,idx),1) + down(z(:,:,idx),1) ...
			     + left(z(:,:,idx),1) + right(z(:,:,idx),1) ...
			   );
		rho(idx) = (rhs(:) \ flat(z(:,:,idx)))/4;
	end
end


