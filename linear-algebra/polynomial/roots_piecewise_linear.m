% Thu 16 Apr 21:29:06 +08 2020
function r = roots_piecewise_linear(cx,cy)
	if (isvector(cx))
		cx = rvec(cx);
	end
	if (isvector(cy))
		cy = rvec(cy);
	end
	if (isvector(cx)&~isvector(cy))
		cx = repmat(cx,size(cy,1),1);
	end
	% linear polynomials
	dy_dx = diff(cy,[],2)./diff(cx,[],2);
	dy_dx(~isfinite(dy_dx)) = 0;
	y0 = cy(:,1:end-1)-cx(:,1:end-1).*dy_dx;
	%dy_dx(abs(dy_dx)<sqrt(eps)) = sqrt(eps);
	r     = NaN(size(cy,1),size(cy,2)-1);
	% check for change of sign
	fdx = (cy(:,1:end-1).*cy(:,2:end) < 0);
	fdx = find(fdx);
	
	% zero
	r(fdx) = -y0(fdx)./dy_dx(fdx);
	% left coincides with zero
	fdx = (0 == cy(:,1:end-1));
	r(fdx) = cx(fdx);
	% last right coincides with zero
	% TODO
	% left and right coincides with zero
	% undefined, set to interval midpoint
%	fdx = ((0 == cy(:,1:end-1)) & (0 == cy(:,2:end)));
%	r(fdx,idx) == 0.5*(cx(fdx);
end

